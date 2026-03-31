% Dan Sadeh 314631763
function [dHbR, dHbO, fig] = CalcNIRS(dataFile, sourceDetectorSep, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
    
    % -------------------------------------------------------------------------
    % Description : The function calculates and plots HbR and HbO concentration 
    %               changes using the Modified Beer-Lambert Law.
    % Assumptions : The data matrix 'd' contains exactly 40 columns (20 channels 
    %               per wavelength). sourceDetectorSep is positive.
    % Inputs :
    %   dataFile - .mat file with intensity data (d), time (t), and Lambda.
    %   sourceDetectorSep - Source-Detector Separation distance in cm.
    %   tissueType - String matching a row in the DPF text file.
    %   plotChannelIdx - Vector indicating channels to plot.
    %   extinctionCoefficientsFile - .csv file with extinction parameters.
    %   DPFperTissueFile - .txt file with baseline DPF multipliers.
    %   relDPFfile - .csv file for wavelength-specific DPF adjustment.
    % Outputs :
    %   dHbR - HbR concentration change matrix (time x 20).
    %   dHbO - HbO concentration change matrix (time x 20).
    %   fig  - Handle to generated figures. Empty if no plots requested.
    % -------------------------------------------------------------------------

    %% Calculation segment - Input Validation and Constants
    fig = [];
    
    % Defining hardware architecture constants to eliminate magic numbers
    NUM_CHANNELS = 20;
    NUM_COLUMNS = 40;
    MICROMOLAR_CONV = 1e6;
    
    if nargin < 3
        error('You must provide at least dataFile, sourceDetectorSep, and tissueType.');
    end
    
    if ~(ischar(dataFile) || isstring(dataFile)) || ~exist(dataFile, 'file')
        error('dataFile must be a valid text string pointing to an existing file.');
    end
    
    if ~isnumeric(sourceDetectorSep) || ~isscalar(sourceDetectorSep) || sourceDetectorSep <= 0
        error('sourceDetectorSep must be a single positive number.');
    end
    
    if nargin < 4 || isempty(plotChannelIdx)
        plotChannelIdx = []; 
    elseif ~isnumeric(plotChannelIdx) || ~isequal(round(plotChannelIdx), plotChannelIdx) || any(plotChannelIdx < 1) || any(plotChannelIdx > NUM_CHANNELS)
        error('plotChannelIdx must contain only whole numbers between 1 and %d.', NUM_CHANNELS);
    end
    
    if nargin < 5 || isempty(extinctionCoefficientsFile); extinctionCoefficientsFile = 'ExtinctionCoefficientsData.csv'; end
    if nargin < 6 || isempty(DPFperTissueFile); DPFperTissueFile = 'DPFperTissue.txt'; end
    if nargin < 7 || isempty(relDPFfile); relDPFfile = 'RelativeDPFCoefficients.csv'; end

    %% Calculation segment - Data Extraction
    dataStruct = load(dataFile);
    
    if ~isfield(dataStruct, 'd') || ~isfield(dataStruct, 't') || ~isfield(dataStruct, 'SD') || ~isfield(dataStruct.SD, 'Lambda')
        error('Loaded file is missing required fields (d, t, SD, or Lambda).');
    end
    
    if size(dataStruct.d, 2) ~= NUM_COLUMNS
        error('Expected the data matrix to have exactly %d columns.', NUM_COLUMNS);
    end
    
    intensityMat = dataStruct.d;
    timeVec = dataStruct.t;
    lambda1 = dataStruct.SD.Lambda(1); 
    lambda2 = dataStruct.SD.Lambda(2); 
    
    extCoeffTable = readtable(extinctionCoefficientsFile);
    dpfTissueTable = readtable(DPFperTissueFile);
    relDpfTable = readtable(relDPFfile);
    
    %% Calculation segment - MBLL Pre-calculations
    % Extracting the tissue-specific pathlength multiplier
    tissueIdx = find(strcmp(dpfTissueTable.Tissue, tissueType));
    if isempty(tissueIdx)
        error('Tissue type not found in DPF table.');
    end
    baseDPF = dpfTissueTable.DPF(tissueIdx);
    
    % Adjusting pathlength based on photon scattering variance at different wavelengths
    lambda1Idx = find(relDpfTable.wavelength == lambda1);
    lambda2Idx = find(relDpfTable.wavelength == lambda2);
    relDPF1 = relDpfTable.relDPFcoeff(lambda1Idx);
    relDPF2 = relDpfTable.relDPFcoeff(lambda2Idx);
    
    % Calculating Effective Pathlength (Leff)
    leffL1 = sourceDetectorSep * baseDPF * relDPF1;
    leffL2 = sourceDetectorSep * baseDPF * relDPF2;
    
    % Defining the inverse (1/Leff) here optimizes loop performance later
    invLeffL1 = 1 / leffL1;
    invLeffL2 = 1 / leffL2;
    
    % Building and immediately inverting the Extinction Coefficient Matrix (Epsilon)
    epsIdxL1 = find(extCoeffTable.wavelength == lambda1);
    epsIdxL2 = find(extCoeffTable.wavelength == lambda2);
    
    epsMatrix = [extCoeffTable.HbO2(epsIdxL1), extCoeffTable.HHb(epsIdxL1); 
                 extCoeffTable.HbO2(epsIdxL2), extCoeffTable.HHb(epsIdxL2)];
    invEpsMatrix = inv(epsMatrix);

    %% Calculation segment - Core MBLL Matrix Multiplication
    % Preallocating arrays limits dynamic memory expansion inside the loop, 
    % significantly reducing computational overhead.
    dHbO = zeros(length(timeVec), NUM_CHANNELS);
    dHbR = zeros(length(timeVec), NUM_CHANNELS);
    
    for channelIdx = 1:NUM_CHANNELS
        % Isolating baseline and continuous data for the specific channel
        intensityZeroL1 = intensityMat(1, channelIdx);
        intensityTimeL1 = intensityMat(:, channelIdx);
        
        intensityZeroL2 = intensityMat(1, NUM_CHANNELS + channelIdx);
        intensityTimeL2 = intensityMat(:, NUM_CHANNELS + channelIdx);
        
        % Calculating Optical Density using base 10 logarithmic attenuation
        opticalDensityL1 = log10(intensityZeroL1 ./ intensityTimeL1);
        opticalDensityL2 = log10(intensityZeroL2 ./ intensityTimeL2);
        
        % Normalizing optical density by the photon's travel distance
        scaledL1 = opticalDensityL1 .* invLeffL1;
        scaledL2 = opticalDensityL2 .* invLeffL2;
        
        % Solving the linear system to extract distinct HbO/HbR chromophore concentrations
        rightSideMatrix = [scaledL1, scaledL2]';
        concentrations = invEpsMatrix * rightSideMatrix; 
        
        dHbO(:, channelIdx) = concentrations(1, :)';
        dHbR(:, channelIdx) = concentrations(2, :)';
    end

    %% Plot segment - Data Visualization
    if isempty(plotChannelIdx)
        disp('No channels were given so no plot is generated.');
    else
        for ii = 1:length(plotChannelIdx)
            channelToPlot = plotChannelIdx(ii); 
            fig(ii) = figure;
            
            % Multiplying by the conversion constant brings the values into a readable microMolar integer range
            plot(timeVec, dHbO(:, channelToPlot) * MICROMOLAR_CONV, 'r', 'LineWidth', 1.5);
            hold on;
            plot(timeVec, dHbR(:, channelToPlot) * MICROMOLAR_CONV, 'b', 'LineWidth', 1.5);
            
            xlim([0, max(timeVec)]);
            title(sprintf('Concentration changes as func of time for channel %d', channelToPlot));
            xlabel('Time (seconds)');
            ylabel('\Delta Concentration (\muM)');
            legend('\Delta HbO', '\Delta HbR', 'Location', 'best');
            grid on;
            hold off;
        end
    end
end