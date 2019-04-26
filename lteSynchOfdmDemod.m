clear;
clc;
% load('Capture-2019.04.02.04.49.47.363.mat')
% load('Capture-2019.04.02.04.49.02.081.mat')
load('Capture-2019.04.02.04.38.27.494.mat');

B = 20e6;
[yf dt] = FilRs(Y,XDelta,B);

sr = 30.72e6;
    
B=20e6;
Bsync = 1.25e6;
[ysync dt_sync] = FilDs(yf,dt,B,Bsync);
      
%==================  HOUSE KEEPING  ===========================================

% Set up some housekeeping variables:
% separator for command window logging
% separator = repmat('-',1,50);

    channelFigure = figure('Visible','off');

[spectrumAnalyzer,synchCorrPlot,pdcchConstDiagram] = ...
    hSIB1RecoveryExamplePlots(channelFigure,sr);
% PDSCH EVM
pdschEVM = comm.EVM();
pdschEVM.MaximumEVMOutputPort = true;

% The sampling rate for the initial cell search is established using
% lteOFDMInfo configured for 6 resource blocks. enb.CyclicPrefix is set
% temporarily in the call to lteOFDMInfo to suppress a default value
% warning (it does not affect the sampling rate).
enb = struct;                   % eNodeB config structure
enb.NDLRB = 6;                  % Number of resource blocks
ofdmInfo = lteOFDMInfo(setfield(enb,'CyclicPrefix','Normal')); %#ok<SFLD>

%==================  CELL SEARCH ==================================================


fprintf('\nPerforming cell search...\n');


    cyclicPrefixes = {'Normal' 'Extended'};

% Perform cell search across cyclic prefix length
% option and record the value with the maximum correlation; if
% multiple cell search is configured, this example will decode the first
% (strongest) detected cell

searchalg.MaxCellCount = 1;
searchalg.SSSDetection = 'PostFFT';
peakMax = -Inf;

 enb.DuplexMode = 'FDD';

    for cyclicPrefix = cyclicPrefixes
        enb.CyclicPrefix = cyclicPrefix{1};
        [enb.NCellID, offset, peak] = lteCellSearch(enb, ysync, searchalg);
        enb.NCellID = enb.NCellID(1);
        offset = offset(1);
        peak = peak(1);
        if (peak>peakMax)
            enbMax = enb;
            offsetMax = offset;
            peakMax = peak;
        end
    end



% Use the cell identity, cyclic prefix length, duplex mode and timing
% offset which gave the maximum correlation during cell search
enb = enbMax;
offset = offsetMax;

% Compute the correlation for each of the three possible primary cell
% identities; the peak of the correlation for the cell identity established
% above is compared with the peak of the correlation for the other two
% primary cell identities in order to establish the quality of the
% correlation.

corr = cell(1,3);
idGroup = floor(enbMax.NCellID/3);
for i = 0:2
    enb.NCellID = idGroup*3 + mod(enbMax.NCellID + i,3);
    [~,corr{i+1}] = lteDLFrameOffset(enb, ysync);
    corr{i+1} = sum(corr{i+1},2);
end
threshold = 1.3 * max([corr{2}; corr{3}]); % multiplier of 1.3 empirically obtained
if (max(corr{1})<threshold)
    warning('Synchronization signal correlation was weak; detected cell identity may be incorrect.');
end
% Return to originally detected cell identity
enb.NCellID = enbMax.NCellID;

% Plot PSS/SSS correlation and threshold
synchCorrPlot.YLimits = [0 max([corr{1}; threshold])*1.1];
synchCorrPlot([corr{1} threshold*ones(size(corr{1}))]);

% Perform timing synchronization
fprintf('Timing offset to frame start: %d samples\n',offset);
ysync = ysync(1+offset:end,:);
enb.NSubframe = 0;

% Show cell-wide settings
fprintf('Cell-wide settings after cell search:\n');
disp(enb);

%===================  Frequency Offset Estimation and Correction =======================

fprintf('\nPerforming frequency offset estimation...\n');

delta_f = lteFrequencyOffset(enb, ysync);
fprintf('Frequency offset: %0.3fHz\n',delta_f);
ysync = lteFrequencyCorrect(enb, ysync, delta_f);

%=========================== OFDM Demodulation and Channel Estimation ================

% Channel estimator configuration
cec.PilotAverage = 'UserDefined';     % Type of pilot averaging
cec.FreqWindow = 13;                  % Frequency window size
cec.TimeWindow = 9;                   % Time window size
cec.InterpType = 'cubic';             % 2D interpolation type
cec.InterpWindow = 'Centered';        % Interpolation window type
cec.InterpWinSize = 1;                % Interpolation window size

% Assume 4 cell-specific reference signals for initial decoding attempt;
% ensures channel estimates are available for all cell-specific reference
% signals
enb.CellRefP = 4;

fprintf('Performing OFDM demodulation...\n\n');

griddims = lteResourceGridSize(enb); % Resource grid dimensions
L = griddims(2);                     % Number of OFDM symbols in a subframe
% OFDM demodulate signal
rxgrid = lteOFDMDemodulate(enb, ysync);
if (isempty(rxgrid))
    fprintf('After timing synchronization, signal is shorter than one subframe so no further demodulation will be performed.\n');
    return;
end
% Perform channel estimation
[hest, nest] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:L,:));

%======================  PBCH Demodulation, BCH Decoding, MIB Parsing ================
% Decode the MIB
% Extract resource elements (REs) corresponding to the PBCH from the first
% subframe across all receive antennas and channel estimates
fprintf('Performing MIB decoding...\n');
pbchIndices = ltePBCHIndices(enb);
[pbchRx, pbchHest] = lteExtractResources( ...
    pbchIndices, rxgrid(:,1:L,:), hest(:,1:L,:,:));

% Decode PBCH
[bchBits, pbchSymbols, nfmod4, mib, enb.CellRefP] = ltePBCHDecode( ...
    enb, pbchRx, pbchHest, nest);

% Parse MIB bits
enb = lteMIB(mib, enb);

% Incorporate the nfmod4 value output from the function ltePBCHDecode, as
% the NFrame value established from the MIB is the System Frame Number
% (SFN) modulo 4 (it is stored in the MIB as floor(SFN/4))
enb.NFrame = enb.NFrame+nfmod4;

% Display cell wide settings after MIB decoding
fprintf('Cell-wide settings after MIB decoding:\n');
disp(enb);

%=============================== OFDM Demodulation on Full Bandwidth ================
fprintf('Restarting reception now that bandwidth (NDLRB=%d) is known...\n',enb.NDLRB);

% Use the full LTE signal bandwidth

waveform_full = yf;

% Perform frequency offset estimation and correction
fprintf('\nPerforming frequency offset estimation...\n');
delta_f = lteFrequencyOffset(enb, waveform_full);
fprintf('Frequency offset: %0.3fHz\n',delta_f);
waveform_full = lteFrequencyCorrect(enb, waveform_full, delta_f);

% Find beginning of frame
fprintf('\nPerforming timing offset estimation...\n');
offset = lteDLFrameOffset(enb, waveform_full);
fprintf('Timing offset to frame start: %d samples\n',offset);
% Aligning signal with the start of the frame
waveform_full = waveform_full(1+offset:end,:);

% OFDM demodulation
fprintf('\nPerforming OFDM demodulation...\n\n');

% This is the array that conatains the received resource elements.
rxgrid = lteOFDMDemodulate(enb, waveform_full);






function [y dt] = FilRs(xt,dti,B)
    N = length(xt);
    df = 1/(N*dti);
    Nb = floor(B/(2*df));
    Fil(1:N) = 0;
    Fil(1:Nb)= 1;
    Fil(N-Nb+1:N) = 1;
    Xf = fft(xt);
    Xff = Fil'.*Xf;
    dt = dti*56/30.72;
    y = resample(double(ifft(Xff)),3072,5600);
end

function [y, dt] = FilDs(xt,dti,Bi,Bo)
    N = length(xt);
    df = 1/(N*dti);
    Nb = floor(Bo/df);
    Fil(1:N) = 0;
    Fil(1:Nb)= 1;
    Fil(N-Nb+1:N) = 1;
    Xf = fft(xt);
    Xff = Fil'.*Xf;
    nds = floor(.1 + Bi/Bo);
    dt = dti*nds;
    y = downsample(double(ifft(Xff)),nds);
end



