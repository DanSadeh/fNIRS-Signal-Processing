% Dan Sadeh 314631763
clear;
clc;
close all;

%% Calculation segment - Initialization and Constants
% Defining hardware and physiological constants upfront to avoid magic numbers
SOURCE_DETECTOR_SEP = 3; 
HEARTBEAT_LOW_HZ = 0.8;
HEARTBEAT_HIGH_HZ = 1.5;
NOISE_FLOOR_HZ = 2.5;
MICROMOLAR_CONV = 1e6;
SECONDS_PER_MINUTE = 60;

% Defining inputs
dataFile1 = '.\FN_032_V1_Postdose1_Nback.mat';
dataFile2 = '.\FN_031_V2_Postdose2_Nback.mat';
tissueType = 'adult_head';
plotChannelIdx = [1, 2];
extinctionCoefficientsFile = '.\ExtinctionCoefficientsData.csv';
DPFperTissueFile = '.\DPFperTissue.txt';
relDPFfile = '.\RelativeDPFCoefficients.csv';

%% Calculation segment - Time Domain Processing (Question 2)
% Applying the MBLL to both recording sessions to extract relative concentrations
[dHbR_1, dHbO_1, fig_1] = CalcNIRS(dataFile1, SOURCE_DETECTOR_SEP, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile);

[dHbR_2, dHbO_2, fig_2] = CalcNIRS(dataFile2, SOURCE_DETECTOR_SEP, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile);

% Automatically exporting Question 2 Time Domain plots at 600 DPI
exportgraphics(figure(fig_1(1)), 'File1_Channel1_TimeDomain.png', 'Resolution', 600);
exportgraphics(figure(fig_1(2)), 'File1_Channel2_TimeDomain.png', 'Resolution', 600);
exportgraphics(figure(fig_2(1)), 'File2_Channel1_TimeDomain.png', 'Resolution', 600);
exportgraphics(figure(fig_2(2)), 'File2_Channel2_TimeDomain.png', 'Resolution', 600);

%% Calculation segment - Frequency Domain Transformation (Question 3)
% Extracting the time vector to determine the hardware sampling rate for the FFT
dataStruct1 = load(dataFile1);
timeVec = dataStruct1.t; 
samplingFreq = 1 / mean(diff(timeVec)); 
numSamples = length(timeVec); 

% Subtracting the mathematical mean eliminates the 0Hz DC offset
centeredSignalChannel1 = dHbO_1(:, 1) - mean(dHbO_1(:, 1)); 

% Executing the discrete Fourier transform
fftResultHbO = fft(centeredSignalChannel1);

% Normalizing by the total sample count ensures the Y-axis reflects physical amplitude
fftAbsHbO = abs(fftResultHbO / numSamples); 

% Slicing the array to discard negative frequencies
fftFinalHbO = fftAbsHbO(1:floor(numSamples/2) + 1); 

% Multiplying by 2 conserves the physical energy
fftFinalHbO(2:end-1) = 2 * fftFinalHbO(2:end-1); 

% Constructing the physical X-axis
freqAxis = (0:floor(numSamples/2)) * (samplingFreq / numSamples); 

%% Calculation segment - SNR & Heart Rate
% We calculate the peak BEFORE plotting so we can explicitly mark it on the graph
heartBeatMask = (freqAxis >= HEARTBEAT_LOW_HZ & freqAxis <= HEARTBEAT_HIGH_HZ);
[heartBeatStrength, maxIdx] = max(fftFinalHbO(heartBeatMask));

% Extracting exact frequency to convert to Beats Per Minute (BPM)
maskedFreqs = freqAxis(heartBeatMask);
heartBeatFreq = maskedFreqs(maxIdx);
bpm = heartBeatFreq * SECONDS_PER_MINUTE; 

% Defining the baseline hardware noise floor
noiseMask = (freqAxis >= NOISE_FLOOR_HZ);
noiseAvg = mean(fftFinalHbO(noiseMask));

% Final ratio calculation
SNR = (heartBeatStrength / noiseAvg);

%% Plot segment - FFT Spectrum Visualization with Peak Marker
fig_fft = figure;
% Magnifying signal to microMolar to view biological changes clearly
plot(freqAxis, fftFinalHbO * MICROMOLAR_CONV, 'g', 'LineWidth', 1.5); 
hold on;

% Marking the peak with a red star and an exact text label (e.g. 57.6 BPM)
plot(heartBeatFreq, heartBeatStrength * MICROMOLAR_CONV, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
text(heartBeatFreq + 0.1, heartBeatStrength * MICROMOLAR_CONV, sprintf('Peak: %.2f Hz (%.1f BPM)', heartBeatFreq, bpm), 'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold');

xlim([0.2, max(freqAxis)]); % Starting from 0.2Hz to exclude slow artifacts
title('Concentration as function of frequency for channel 1');
xlabel('Frequency (Hz)')
ylabel('\Delta Concentration (\muM)');
legend('\Delta HbO', 'Heartbeat Peak', 'location', 'best');
grid on;
hold off;

% Automatically exporting Question 3 Frequency Domain plot at 600 DPI
exportgraphics(fig_fft, 'FFT_Spectrum_Channel1.png', 'Resolution', 600);

% Printouts
fprintf('Heartbeat Signal Strength: %.4e\n', heartBeatStrength);
fprintf('Baseline Noise Level:      %.4e\n', noiseAvg);
fprintf('Heart rate (BPM):          %.4f\n', bpm);
fprintf('Final SNR:                 %.4f\n', SNR);