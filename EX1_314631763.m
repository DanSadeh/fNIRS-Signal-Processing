% Dan Sadeh 314631763
clear;
clc;
close all;

%% Calculation segment - Initialization and Constants
% Defining hardware and physiological constants
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

%% Calculation segment - Frequency Domain Transformation (Question 3)
% Extracting the time vector to determine the hardware sampling rate for the FFT
dataStruct1 = load(dataFile1);
timeVec = dataStruct1.t; 
samplingFreq = 1 / mean(diff(timeVec)); 
numSamples = length(timeVec); 

% Subtracting the mathematical mean eliminates the 0Hz DC offset, preventing 
% the baseline concentration from ruining the amplitude scale of our graph.
centeredSignalChannel1 = dHbO_1(:, 1) - mean(dHbO_1(:, 1)); 

% Executing the discrete Fourier transform
fftResultHbO = fft(centeredSignalChannel1);

% Normalizing by the total sample count ensures the Y-axis reflects true physical 
% concentration changes rather than arbitrary summation energy.
fftAbsHbO = abs(fftResultHbO / numSamples); 

% Slicing the array to discard negative frequencies, as they are mathematical 
% mirrors and hold no biological relevance.
fftFinalHbO = fftAbsHbO(1:floor(numSamples/2) + 1); 

% Multiplying by 2 conserves the physical energy lost from discarding the negative frequency bins.
fftFinalHbO(2:end-1) = 2 * fftFinalHbO(2:end-1); 

% Constructing the physical X-axis by multiplying bin indices by the frequency resolution
freqAxis = (0:floor(numSamples/2)) * (samplingFreq / numSamples); 

%% Plot 1 segment - FFT Spectrum Visualization
figure;
% Magnifying signal to microMolar to view biological changes clearly
plot(freqAxis, fftFinalHbO * MICROMOLAR_CONV, 'g', 'LineWidth', 1.5); 
xlim([0.2, max(freqAxis)]); % Starting from 0.2Hz to exclude huge Mayer wave/respiration noise
title('Concentration as function of frequency for channel 1');
xlabel('Frequency (Hz)')
ylabel('\Delta Concentration (\muM)');
legend('\Delta HbO', 'location', 'best');
grid on;

%% Calculation segment - SNR & Heart Rate
% Isolating the frequency bins that correspond exclusively to the biological cardiac cycle
heartBeatMask = (freqAxis >= HEARTBEAT_LOW_HZ & freqAxis <= HEARTBEAT_HIGH_HZ);
[heartBeatStrength, maxIdx] = max(fftFinalHbO(heartBeatMask));

% Defining the baseline hardware noise floor
noiseMask = (freqAxis >= NOISE_FLOOR_HZ);
noiseAvg = mean(fftFinalHbO(noiseMask));

% Final ratio calculation
SNR = (heartBeatStrength / noiseAvg);

% Extracting exact frequency to convert to Beats Per Minute (BPM)
maskedFreqs = freqAxis(heartBeatMask);
heartBeatFreq = maskedFreqs(maxIdx);
bpm = heartBeatFreq * SECONDS_PER_MINUTE; 

fprintf('Heartbeat Signal Strength: %.4e\n', heartBeatStrength);
fprintf('Baseline Noise Level:      %.4e\n', noiseAvg);
fprintf('Heart rate (BPM):          %.4f\n', bpm);
fprintf('Final SNR:                 %.4f\n', SNR);