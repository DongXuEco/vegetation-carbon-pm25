% Calculate pixel-wise OLS slopes with an intercept using actual years.
% Inputs must share the grid and units; only complete time series are used.
% Set noDataValues for numeric missing-data codes.
% trendPeriodYears = 10 gives change per decade; 1 gives change per year.
clc;
clear;

% 1. Input and output settings
inputFolder = '';       % Folder containing the time-series GeoTIFF files.
filePrefix = '';        % Filename prefix preceding the year.
outputFile = '';        % Full path of the output GeoTIFF file.
years = 2025:5:2100;     % Use 2001:2022 for annual historical data.
trendPeriodYears = 10;   % Set to 1 for annual change or 10 for decadal change.
noDataValues = [];      % Numeric missing-data codes, e.g. [-9999, -32768].

if isempty(inputFolder) || isempty(outputFile)
    error('Set inputFolder and outputFile before running the script.');
end
if numel(years) < 2 || any(~isfinite(years)) || any(diff(years) <= 0)
    error('years must contain at least two strictly increasing finite years.');
end
if ~isscalar(trendPeriodYears) || ~isfinite(trendPeriodYears) || trendPeriodYears <= 0
    error('trendPeriodYears must be a positive finite scalar.');
end

% 2. Read the reference raster and time series
referenceFile = fullfile(inputFolder, [filePrefix, int2str(years(1)), '.tif']);
[a, R] = geotiffread(referenceFile);
info = geotiffinfo(referenceFile);
if ~ismatrix(a)
    error('Input rasters must contain a single band.');
end
[m, n] = size(a);
numYears = numel(years);
datasum = nan(m * n, numYears);

k = 1;
for year = years
    filename = fullfile(inputFolder, [filePrefix, int2str(year), '.tif']);
    [data, currentR] = geotiffread(filename);
    currentInfo = geotiffinfo(filename);
    if ~isequal(size(data), [m, n]) || ~isequal(currentR, R) || ...
            ~isequaln(currentInfo.GeoTIFFCodes, info.GeoTIFFCodes)
        error('Input rasters must use the same grid and coordinate system: %s', filename);
    end
    data = double(data);
    data(~isfinite(data) | ismember(data, noDataValues)) = NaN;
    datasum(:, k) = reshape(data, m * n, 1);
    k = k + 1;
end

% 3. Calculate the OLS slope at each pixel
trend = nan(m, n);
centeredYears = double(years(:)') - mean(double(years(:)));
sumSquaredYears = sum(centeredYears .^ 2);
for i = 1:size(datasum, 1)
    data = datasum(i, :);
    if all(isfinite(data))
        centeredData = data - mean(data);
        trend(i) = sum(centeredYears .* centeredData) / sumSquaredYears ...
            * trendPeriodYears;
    end
end

% 4. Export the trend raster with the reference spatial metadata
geotiffwrite(outputFile, trend, R, ...
    'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
fprintf('Trend raster saved to: %s\n', outputFile);
