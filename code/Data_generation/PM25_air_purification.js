/*
 * Annual vegetation PM2.5 adsorption: Q = P * Vd * LAI * T / 1e6.
 * P: micrograms/m3; Vd: m/s; LAI: m2/m2; Q: g/m2/year.
 * T = days with precipitation below 0.1 mm/day * 86400 seconds.
 * Use the Vd band for the corresponding year.
 */
// Common valid vegetation domain for all variables, years and scenarios.
// Single-band raster: 1 = included, 0 = excluded; use the same asset in all scripts.
var vegetationMaskAsset = '';

if (!vegetationMaskAsset) {
    throw new Error('Fill in the common vegetation mask asset ID.');
}


// Annual PM2.5 asset IDs (single band, micrograms/m3).
var pm25Assets = {
    2001: ''
};
// Multiyear deposition velocity raster asset (m/s), with one band per year.
// Band names: Vd_2001, Vd_2002, ...; adjust vdBandPrefix to match the asset.
var vdAsset = '';
var vdBandPrefix = 'Vd_';

var startYear = 2001;
var endYear = 2001;
var scale = 10000;
var exportResults = false;

// Study area
var geometry = 
    ee.Geometry.MultiPolygon(
        [[[[-179.31132812500005, 82.57599102151069],
          [-179.31132812500005, -67.37002275111178],
          [-168.588671875, -67.37002275111178],
          [-168.588671875, 82.57599102151069]]],
        [[[-170.17070312500002, 82.54184417673481],
          [-170.17070312500002, -67.4712617176781],
          [181.74335937499995, -67.4712617176781],
          [181.74335937499995, 82.54184417673481]]]], null, false);

var mask = ee.Image(vegetationMaskAsset).eq(1).unmask(0)
    .clip(geometry).setDefaultProjection('EPSG:4326', null, scale);

if (!vdAsset) {
    throw new Error('Missing deposition velocity asset ID.');
}
var vdData = ee.Image(vdAsset);

for (var i = startYear; i <= endYear; i++) {
    var yearOut = String(i);
    var yearOut2 = String(i + 1);
    if (!pm25Assets[i]) {
        throw new Error('Missing PM2.5 asset ID for ' + yearOut);
    }

    // Deposition velocity for the current year.
    var V_pm25 = vdData.select(vdBandPrefix + yearOut)
        .clip(geometry).toFloat().setDefaultProjection('EPSG:4326', null, scale);

    // LAI
    var LAI = ee.ImageCollection('MODIS/061/MOD15A2H')
        .filterDate(yearOut + '-01-01', yearOut2 + '-01-01')
        .filterBounds(geometry).select('Lai_500m')
        .mean().divide(10)
        .clip(geometry).unmask(0).toFloat()
        .setDefaultProjection('EPSG:4326', null, scale);

    // Time below the precipitation threshold (s/year)
    var pre = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filterDate(yearOut + '-01-01', yearOut2 + '-01-01')
        .select('total_precipitation_sum');
    var preThreshold = 0.0001; // m/day = 0.1 mm/day
    var preBelowThreshold = pre.map(function(image) {
        var preMask = image.lt(preThreshold).unmask(0);
        return preMask.rename('pre_below_threshold')
            .copyProperties(image, ['system:time_start']);
    });
    var pretimeCount = preBelowThreshold.sum().multiply(86400.0)
        .clip(geometry).toFloat().setDefaultProjection('EPSG:4326', null, scale);

    var PM25 = ee.Image(pm25Assets[i]).clip(geometry).unmask(0)
        .multiply(mask).toFloat().setDefaultProjection('EPSG:4326', null, scale);
    var PM25_NEW = PM25.updateMask(mask);

    // PM2.5 adsorption (g/m2/year)
    var a = V_pm25.multiply(LAI)
        .multiply(pretimeCount).clip(geometry).toFloat()
        .setDefaultProjection('EPSG:4326', null, scale);
    var Q_pm25 = PM25_NEW.multiply(a).divide(1000000.0).updateMask(mask);

    var vis = {min: 0, max: 30,
        palette: ['ffffff', 'b7f0ae', '21f600', '0000FF', 'FDFF92', 'FF2700', 'd600ff']};
    Map.addLayer(Q_pm25, vis, 'PM2.5 adsorption ' + yearOut);

    if (exportResults) {
        Export.image.toDrive({
            image: Q_pm25,
            description: 'PM25_adsorption_' + yearOut,
            folder: 'PM25_adsorption',
            region: geometry,
            crs: 'EPSG:4326',
            scale: scale,
            maxPixels: 1e13
        });
    }
}
