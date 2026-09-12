/*
 * SSP adsorption: Q = P * Vd_2022 * LAI * T_2022 / 1e6.
 * P: micrograms/m3; Vd: m/s; T: seconds/year; LAI: physical, unscaled values.
 * Vd and T remain fixed at 2022; P and LAI vary by scenario and period.
 * Output is g/m2/year, not a five-year total.
 */
// Common valid vegetation domain for all variables, years and scenarios.
// Single-band raster: 1 = included, 0 = excluded; use the same asset in all scripts.
var vegetationMaskAsset = '';

if (!vegetationMaskAsset) {
    throw new Error('Fill in the common vegetation mask asset ID.');
}


var vdAsset = '';       // Single-band deposition velocity raster for 2022 (m/s).
var noPreTimeAsset = ''; // Single-band no-precipitation duration for 2022 (s/year).

// PM2.5 concentration (micrograms/m3).
var pm25Assets = {
    SSP126: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    },
    SSP245: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    },
    SSP585: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    }
};

// LAI in physical units (m2/m2).
var laiAssets = {
    SSP126: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    },
    SSP245: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    },
    SSP585: {
        2025: '',
        2030: '',
        2035: '',
        2040: '',
        2045: '',
        2050: '',
        2055: '',
        2060: '',
        2065: '',
        2070: '',
        2075: '',
        2080: '',
        2085: '',
        2090: '',
        2095: '',
        2100: ''
    }
};

var scenarios = ['SSP126', 'SSP245', 'SSP585'];
var startYear = 2025;
var endYear = 2100;
var yearStep = 5;
var scale = 10000;
var exportResults = false;

if (!vdAsset || !noPreTimeAsset) {
    throw new Error('Fill in the 2022 Vd and no-precipitation duration asset IDs.');
}
for (var s = 0; s < scenarios.length; s++) {
    var scenario = scenarios[s];
    for (var y = startYear; y <= endYear; y += yearStep) {
        if (!pm25Assets[scenario] || !pm25Assets[scenario][y] ||
            !laiAssets[scenario] || !laiAssets[scenario][y]) {
            throw new Error('Missing PM2.5 or LAI asset for ' + scenario + ' ' + y);
        }
    }
}

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

// Fixed 2022 deposition background.
var V_pm25 = ee.Image(vdAsset).clip(geometry).toFloat()
    .setDefaultProjection('EPSG:4326', null, scale);
var pretimeCount = ee.Image(noPreTimeAsset).clip(geometry).toFloat()
    .setDefaultProjection('EPSG:4326', null, scale);

for (var s = 0; s < scenarios.length; s++) {
    var scenario = scenarios[s];
    for (var i = startYear; i <= endYear; i += yearStep) {
        var yearOut = String(i);

        // LAI
        var LAI = ee.Image(laiAssets[scenario][i])
            .clip(geometry).unmask(0).toFloat()
            .setDefaultProjection('EPSG:4326', null, scale);

        var PM25 = ee.Image(pm25Assets[scenario][i]).clip(geometry).unmask(0)
            .multiply(mask).toFloat().setDefaultProjection('EPSG:4326', null, scale);
        var PM25_NEW = PM25.updateMask(mask);

        // PM2.5 adsorption (g/m2/year)
        var a = V_pm25.multiply(LAI)
            .multiply(pretimeCount).clip(geometry).toFloat()
            .setDefaultProjection('EPSG:4326', null, scale);
        var Q_pm25 = PM25_NEW.multiply(a).divide(1000000.0).updateMask(mask);

        var vis = {min: 0, max: 30,
            palette: ['ffffff', 'b7f0ae', '21f600', '0000FF', 'FDFF92', 'FF2700', 'd600ff']};
        Map.addLayer(Q_pm25, vis, 'PM2.5 adsorption ' + scenario + ' ' + yearOut);

        if (exportResults) {
            Export.image.toDrive({
                image: Q_pm25,
                description: 'PM25_adsorption_' + scenario + '_' + yearOut,
                folder: 'PM25_adsorption',
                region: geometry,
                crs: 'EPSG:4326',
                scale: scale,
                maxPixels: 1e13
            });
        }
    }
}
