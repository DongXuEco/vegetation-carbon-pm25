/*
 * Local Pearson correlation and two-sided P-values in 19 x 19 pixel windows.
 * Radius: 9 pixels on a nominal 10 km grid; ground distances depend on projection.
 * Use carbon and adsorption rasters from the same period and vegetation domain.
 * Windows need at least three paired values and nonzero variance.
 * P-values are not corrected for spatial autocorrelation or multiple testing.
 */

// 1. Study region and input assets
var roi = ee.Geometry.MultiPolygon(
    [[[[-179.31132812500005, 82.57599102151069],
       [-179.31132812500005, -67.37002275111178],
       [-168.588671875, -67.37002275111178],
       [-168.588671875, 82.57599102151069]]],
     [[[-170.17070312500002, 82.54184417673481],
       [-170.17070312500002, -67.4712617176781],
       [181.74335937499995, -67.4712617176781],
       [181.74335937499995, 82.54184417673481]]]], null, false);

var pm25Asset = '';   // PM2.5 adsorption raster for the selected period/scenario.
var carbonAsset = ''; // Carbon stock raster for the same period/scenario.
var exportFolder = ''; // Google Drive output folder.
var analysisScale = 10000;

if (!pm25Asset || !carbonAsset || !exportFolder) {
  throw new Error('Set the input asset IDs and the Google Drive output folder.');
}

var pm25 = ee.Image(pm25Asset).select([0]).rename('pm25').toDouble();
var carbon = ee.Image(carbonAsset).select([0]).rename('carbon').toDouble();

// 2. Align the inputs and retain paired valid pixels
var analysisProjection = carbon.projection().atScale(analysisScale);
pm25 = pm25.reproject({crs: analysisProjection});
carbon = carbon.reproject({crs: analysisProjection});

var commonMask = pm25.mask().gt(0).and(carbon.mask().gt(0));
pm25 = pm25.updateMask(commonMask);
carbon = carbon.updateMask(commonMask);

// Use the same grid for neighborhood calculations and exports.
var exportGrid = analysisProjection.getInfo();

// 3. Calculate local Pearson correlation and its two-sided p-value
var computeCorrelation = function(image1, image2, kernel) {
  var paired = image1.rename('pm25').addBands(image2.rename('carbon'));

  var count = paired.select('pm25').reduceNeighborhood({
    reducer: ee.Reducer.count(),
    kernel: kernel,
    skipMasked: true
  });

  var ranges = paired.reduceNeighborhood({
    reducer: ee.Reducer.minMax(),
    kernel: kernel,
    skipMasked: true
  });

  var validWindow = count.gte(3)
      .and(ranges.select('pm25_max').gt(ranges.select('pm25_min')))
      .and(ranges.select('carbon_max').gt(ranges.select('carbon_min')));

  var result = paired.reduceNeighborhood({
    reducer: ee.Reducer.pearsonsCorrelation()
        .setOutputs(['correlation', 'p_value']),
    kernel: kernel,
    skipMasked: true
  });

  return result.updateMask(validWindow);
};

for (var windowSize = 19; windowSize < 20; windowSize = windowSize + 2) {
  var flag = String(windowSize);
  var kernel = ee.Kernel.square({
    radius: (windowSize - 1) / 2,
    units: 'pixels',
    normalize: false
  });

  var result = computeCorrelation(pm25, carbon, kernel);
  var correlationImage = result.select('correlation');
  var pValueImage = result.select('p_value');

  // 4. Export correlation and p-value rasters
  Export.image.toDrive({
    folder: exportFolder,
    image: correlationImage,
    description: 'CorrelationImage' + flag,
    fileNamePrefix: 'CorrelationImage' + flag,
    crs: exportGrid.crs,
    crsTransform: exportGrid.transform,
    region: roi,
    maxPixels: 1e15,
    fileFormat: 'GeoTIFF'
  });

  Export.image.toDrive({
    folder: exportFolder,
    image: pValueImage,
    description: 'PValueImage' + flag,
    fileNamePrefix: 'PValueImage' + flag,
    crs: exportGrid.crs,
    crsTransform: exportGrid.transform,
    region: roi,
    maxPixels: 1e15,
    fileFormat: 'GeoTIFF'
  });
}
