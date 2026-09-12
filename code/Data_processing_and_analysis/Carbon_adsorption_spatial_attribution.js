/*
 * Identify LAI or PM2.5 drivers of local Q-C spatial correlation.
 * Use period means from the same scenario, in physical units and a common mask.
 * 11/12: LAI dominant, positive/negative R; 21/22: PM2.5 dominant, positive/negative R.
 * First-order approximation: omits interaction and assumes positive, locally constant Vd*T.
 * Nominal 10 km geographic pixels vary in ground size with latitude.
 * P-values are uncorrected for spatial autocorrelation and multiple testing.
 */

// 1. Study region and input settings
var roi = ee.Geometry.MultiPolygon(
    [[[[-179.31132812500005, 82.57599102151069],
       [-179.31132812500005, -67.37002275111178],
       [-168.588671875, -67.37002275111178],
       [-168.588671875, 82.57599102151069]]],
     [[[-170.17070312500002, 82.54184417673481],
       [-170.17070312500002, -67.4712617176781],
       [181.74335937499995, -67.4712617176781],
       [181.74335937499995, 82.54184417673481]]]], null, false);

var carbonAsset = '';
var adsorptionAsset = '';
var laiAsset = '';
var pm25Asset = '';
var exportFolder = '';
var scenario = '585';
var windowSize = 9;
var analysisScale = 10000;
var alpha = 0.05;
var applySignificanceMask = true; // false retains all otherwise valid classes.

if (!carbonAsset || !adsorptionAsset || !laiAsset || !pm25Asset || !exportFolder) {
  throw new Error('Set the four input asset IDs and the Google Drive output folder.');
}
if (windowSize < 3 || windowSize % 2 !== 1) {
  throw new Error('windowSize must be an odd integer of at least 3.');
}
if (!(analysisScale > 0) || !(alpha > 0 && alpha < 1)) {
  throw new Error('Use a positive analysisScale and alpha between 0 and 1.');
}

// 2. Align inputs and apply their common valid vegetation mask
var analysisProjection = ee.Projection('EPSG:4326').atScale(analysisScale);
var readInput = function(asset, bandName) {
  var image = ee.Image(asset).select([0]).rename(bandName).toDouble();
  return image.updateMask(image.gte(0)).reproject({crs: analysisProjection}).clip(roi);
};
var Carbon = readInput(carbonAsset, 'carbon');
var Q_absorption = readInput(adsorptionAsset, 'adsorption');
var LAI = readInput(laiAsset, 'lai');
var PM25_conc = readInput(pm25Asset, 'pm25');
var inputs = Carbon.addBands(Q_absorption).addBands(LAI).addBands(PM25_conc);
var commonMask = inputs.mask().reduce(ee.Reducer.min()).gt(0);
inputs = inputs.updateMask(commonMask);
Carbon = inputs.select('carbon');
Q_absorption = inputs.select('adsorption');
LAI = inputs.select('lai');
PM25_conc = inputs.select('pm25');

var kernel = ee.Kernel.square({
  radius: (windowSize - 1) / 2,
  units: 'pixels',
  normalize: false
});
var neighborhoodMean = function(image) {
  return image.reduceNeighborhood({
    reducer: ee.Reducer.mean(), kernel: kernel, skipMasked: true
  });
};

// 3. Calculate covariance contributions on the same neighborhood samples
var mean_C = neighborhoodMean(Carbon);
var mean_LAI = neighborhoodMean(LAI);
var mean_PM25 = neighborhoodMean(PM25_conc);
var cov_LAIC = neighborhoodMean(LAI.multiply(Carbon))
    .subtract(mean_LAI.multiply(mean_C));
var cov_PMC = neighborhoodMean(PM25_conc.multiply(Carbon))
    .subtract(mean_PM25.multiply(mean_C));
var con_LAI = mean_PM25.multiply(cov_LAIC);
var con_PM25 = mean_LAI.multiply(cov_PMC);

// 4. Calculate Pearson correlation and its two-sided p-value
var correlation = Q_absorption.addBands(Carbon).reduceNeighborhood({
  reducer: ee.Reducer.pearsonsCorrelation().setOutputs(['correlation_R', 'p_value']),
  kernel: kernel,
  skipMasked: true
});
var sampleCount = Carbon.reduceNeighborhood({
  reducer: ee.Reducer.count(), kernel: kernel, skipMasked: true
});
var ranges = Q_absorption.addBands(Carbon).reduceNeighborhood({
  reducer: ee.Reducer.minMax(), kernel: kernel, skipMasked: true
});
var validWindow = sampleCount.gte(3)
    .and(ranges.select('adsorption_max').gt(ranges.select('adsorption_min')))
    .and(ranges.select('carbon_max').gt(ranges.select('carbon_min')));
var R_map = correlation.select('correlation_R').updateMask(validWindow);
var pValue = correlation.select('p_value').updateMask(validWindow);
var p_mask = pValue.lt(alpha);
var significant_R = R_map.updateMask(p_mask);

// 5. Preserve the four original attribution classes
var abs_con_LAI = con_LAI.abs();
var abs_con_PM25 = con_PM25.abs();
var dominantMap = R_map.multiply(0)
    .where(abs_con_LAI.gt(abs_con_PM25).and(R_map.gt(0)), 11)
    .where(abs_con_LAI.gt(abs_con_PM25).and(R_map.lt(0)), 12)
    .where(abs_con_PM25.gt(abs_con_LAI).and(R_map.gt(0)), 21)
    .where(abs_con_PM25.gt(abs_con_LAI).and(R_map.lt(0)), 22)
    .rename('attribution');
dominantMap = dominantMap.updateMask(dominantMap.neq(0));
if (applySignificanceMask) {
  dominantMap = dominantMap.updateMask(p_mask);
}
dominantMap = dominantMap.toInt16();

// 6. Display correlation and discrete attribution classes
var rVis = {min: -1, max: 1, palette: ['#0000FF', '#FFFFFF', '#FF0000']};
Map.addLayer(R_map, rVis, 'Q-C spatial correlation', false);
Map.addLayer(significant_R, rVis, 'Q-C correlation: nominal P < ' + alpha, false);
var displayClasses = dominantMap.remap([11, 12, 21, 22], [1, 2, 3, 4]);
Map.addLayer(displayClasses, {
  min: 1, max: 4, palette: ['orange', 'blue', 'green', 'yellow']
}, 'Spatial attribution: 11, 12, 21, 22');

// 7. Export both results on the analysis grid
var exportGrid = analysisProjection.getInfo();
var suffix = scenario + '_' + windowSize + 'x' + windowSize;
Export.image.toDrive({
  image: significant_R.toFloat(),
  folder: exportFolder,
  description: 'Spatial_Significant_R_' + suffix,
  fileNamePrefix: 'Spatial_Significant_R_' + suffix,
  crs: exportGrid.crs,
  crsTransform: exportGrid.transform,
  region: roi,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: dominantMap,
  folder: exportFolder,
  description: 'Spatial_Attribution_' + suffix,
  fileNamePrefix: 'Spatial_Attribution_' + suffix,
  crs: exportGrid.crs,
  crsTransform: exportGrid.transform,
  region: roi,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
