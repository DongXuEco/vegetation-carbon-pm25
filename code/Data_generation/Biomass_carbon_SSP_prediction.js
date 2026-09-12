/*
 * Train biomass models once and predict the three SSP scenarios.
 * Future inputs are prepared five-model means for each five-year period.
 * Predictor bands and units must match training; outputs use the target units.
 */
// Common valid vegetation domain for all variables, years and scenarios.
// Single-band raster: 1 = included, 0 = excluded; use the same asset in all scripts.
var vegetationMaskAsset = '';

if (!vegetationMaskAsset) {
    throw new Error('Fill in the common vegetation mask asset ID.');
}


var predictorsAsset = ''; // Reference-year multiband environmental image.
var samplesAsset = '';    // Point samples with biomass carbon density values.
var targetProperty = 'biomass';
var selectedBands = [];   // Empty: use all predictor bands; exclude targets and IDs.
var referenceYear = 2010;
var scale = 10000;
var crs = 'EPSG:4326';
var split = 0.8;
var seed = 2;
var innerFolds = 10;
var exportResults = true;
var exportFolder = 'Biomass_carbon_SSP';

if (!predictorsAsset || !samplesAsset) {
    throw new Error('Fill in predictorsAsset and samplesAsset.');
}

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

var vegetationMask = ee.Image(vegetationMaskAsset).eq(1).unmask(0)
    .clip(geometry).setDefaultProjection(crs, null, scale);

var predictors = ee.Image(predictorsAsset);
var bands = selectedBands.length ? ee.List(selectedBands) : predictors.bandNames();
predictors = predictors.select(bands).clip(geometry).toFloat()
    .setDefaultProjection(crs, null, scale);
var samplePoints = ee.FeatureCollection(samplesAsset).filterBounds(geometry)
    .filter(ee.Filter.notNull([targetProperty]))
    .map(function(feature) { return feature.set('_sample_id', feature.id()); });
var training = predictors.sampleRegions({
    collection: samplePoints,
    properties: [targetProperty, '_sample_id'],
    projection: ee.Projection(crs).atScale(scale),
    scale: scale,
    tileScale: 16,
    geometries: true
}).filter(ee.Filter.notNull(bands.add(targetProperty))).map(function(feature) {
    var xy = feature.geometry().coordinates();
    var pixel = ee.Number(xy.get(0)).format('%.8f').cat('_')
        .cat(ee.Number(xy.get(1)).format('%.8f'));
    return feature.set('_pixel_id', pixel);
});

// Samples in the same analysis pixel stay in the same partition.
var withRandom = training.randomColumn('_split', seed, 'uniform', ['_pixel_id']);
var trainingData = withRandom.filter(ee.Filter.lt('_split', split));
var validationData = withRandom.filter(ee.Filter.gte('_split', split));
print('Predictor bands', bands);
print('Input / usable / development / test samples', samplePoints.size(),
    training.size(), trainingData.size(), validationData.size());

var modelNames = ['RF', 'GBT', 'CART'];
var bands2 = ['predicted_RF', 'predicted_GBT', 'predicted_CART'];
var bands3 = ['predicted_RF_final', 'predicted_GBT_final', 'predicted_CART_final'];
var outputNames = ['primary_RF', 'primary_GBT', 'primary_CART',
    'secondary_RF', 'secondary_GBT', 'secondary_CART'];

function assignFolds(data, count, foldSeed, column) {
    return data.randomColumn(column + '_random', foldSeed, 'uniform', ['_pixel_id'])
        .map(function(feature) {
            return feature.set(column,
                ee.Number(feature.get(column + '_random')).multiply(count).floor());
        });
}

function trainModels(data, inputBands) {
    return {
        RF: ee.Classifier.smileRandomForest(100).setOutputMode('REGRESSION').train({
            features: data, classProperty: targetProperty, inputProperties: inputBands
        }),
        GBT: ee.Classifier.smileGradientTreeBoost({numberOfTrees: 200, shrinkage: 0.1})
            .setOutputMode('REGRESSION').train({
                features: data, classProperty: targetProperty, inputProperties: inputBands
            }),
        CART: ee.Classifier.smileCart().setOutputMode('REGRESSION').train({
            features: data, classProperty: targetProperty, inputProperties: inputBands
        })
    };
}

function predictTable(data, models, names) {
    return data.classify(models.RF, names[0])
        .classify(models.GBT, names[1]).classify(models.CART, names[2]);
}

function predictImage(data, models, names) {
    return data.classify(models.RF, names[0])
        .addBands(data.classify(models.GBT, names[1]))
        .addBands(data.classify(models.CART, names[2]));
}

// First stacking level: train meta-models on base-model out-of-fold predictions.
function fitPrimary(data, fitSeed) {
    var folded = assignFolds(data, innerFolds, fitSeed, '_base_fold');
    var oof = ee.FeatureCollection([]);
    for (var f = 0; f < innerFolds; f++) {
        var train = folded.filter(ee.Filter.neq('_base_fold', f));
        var test = folded.filter(ee.Filter.eq('_base_fold', f));
        var fitted = trainModels(train, bands);
        oof = oof.merge(predictTable(test, fitted, bands2));
    }
    return {base: trainModels(data, bands), primary: trainModels(oof, bands2)};
}

// Second stacking level: refit the entire first level inside each inner fold.
function fitStack(data, fitSeed) {
    var folded = assignFolds(data, innerFolds, fitSeed, '_meta_fold');
    var oof = ee.FeatureCollection([]);
    for (var f = 0; f < innerFolds; f++) {
        var train = folded.filter(ee.Filter.neq('_meta_fold', f));
        var test = folded.filter(ee.Filter.eq('_meta_fold', f));
        var fitted = fitPrimary(train, fitSeed + 101 + f);
        var basePredictions = predictTable(test, fitted.base, bands2);
        oof = oof.merge(predictTable(basePredictions, fitted.primary, bands3));
    }
    var primaryFit = fitPrimary(data, fitSeed + 1001);
    return {base: primaryFit.base, primary: primaryFit.primary,
        secondary: trainModels(oof, bands3)};
}

function predictStackTable(data, fitted) {
    var first = predictTable(data, fitted.base, bands2);
    var second = predictTable(first, fitted.primary, bands3);
    var third = predictTable(second, fitted.secondary, outputNames.slice(3));
    return third.map(function(feature) {
        return feature.set({primary_RF: feature.get(bands3[0]),
            primary_GBT: feature.get(bands3[1]), primary_CART: feature.get(bands3[2])});
    });
}

function predictStackImage(image, fitted, name) {
    var first = predictImage(image.select(bands), fitted.base, bands2);
    var second = predictImage(first, fitted.primary, bands3);
    return second.classify(fitted.secondary[name], 'carbon_stock')
        .updateMask(vegetationMask).clip(geometry).toFloat().setDefaultProjection(crs, null, scale);
}

function metrics(data, prediction, label, fold) {
    var meanObserved = ee.Number(data.aggregate_mean(targetProperty));
    var residuals = data.map(function(feature) {
        var observed = ee.Number(feature.get(targetProperty));
        var residual = ee.Number(feature.get(prediction)).subtract(observed);
        return feature.set({_error: residual, _squared_error: residual.pow(2),
            _total_squared: observed.subtract(meanObserved).pow(2)});
    });
    var sse = ee.Number(residuals.aggregate_sum('_squared_error'));
    var sst = ee.Number(residuals.aggregate_sum('_total_squared'));
    var correlation = data.reduceColumns(ee.Reducer.pearsonsCorrelation(),
        [targetProperty, prediction]).get('correlation');
    var regression = data.reduceColumns(ee.Reducer.linearFit(),
        [targetProperty, prediction]);
    return ee.Feature(null, {
        dataset: label, fold: fold, model: prediction, n: data.size(),
        R2: ee.Algorithms.If(sst.gt(0), ee.Number(1).subtract(sse.divide(sst)), null),
        r_squared: ee.Algorithms.If(ee.Algorithms.IsEqual(correlation, null),
            null, ee.Number(correlation).pow(2)),
        bias: residuals.aggregate_mean('_error'),
        RMSE: ee.Number(residuals.aggregate_mean('_squared_error')).sqrt(),
        slope: regression.get('scale'), intercept: regression.get('offset')
    });
}

function exportImage(image, name) {
    if (exportResults) {
        Export.image.toDrive({image: image, folder: exportFolder,
            description: name, fileNamePrefix: name,
            region: geometry, crs: crs, scale: scale, maxPixels: 1e13});
    }
}

// Fit the historical models before predicting any scenario or year.
var fittedFinal = fitStack(trainingData, seed + 1000000);
var testPredictions = predictStackTable(validationData, fittedFinal);
var testMetrics = ee.FeatureCollection([
    metrics(testPredictions, 'secondary_RF', 'independent_test', -1)
]);
if (exportResults) {
    Export.table.toDrive({collection: testMetrics, folder: exportFolder,
        description: 'Biomass_SSP_test_metrics', fileFormat: 'CSV'});
    Export.table.toDrive({collection: testPredictions, folder: exportFolder,
        description: 'Biomass_SSP_test_predictions', fileFormat: 'CSV',
        selectors: ['_sample_id', '_pixel_id', targetProperty, 'secondary_RF']});
}

function predictBiomassSSP() {
    var trained_RF = fittedFinal.base.RF;
    var trained_GBT = fittedFinal.base.GBT;
    var trained_CART = fittedFinal.base.CART;
    var trained_RF_meta = fittedFinal.primary.RF;
    var trained_GBT_meta = fittedFinal.primary.GBT;
    var trained_CART_meta = fittedFinal.primary.CART;
    var trained_RF_meta_final = fittedFinal.secondary.RF;

    var predictionAssets = {
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

    for (var s = 0; s < scenarios.length; s++) {
        var scenario = scenarios[s];
        for (var y = startYear; y <= endYear; y += yearStep) {
            if (!predictionAssets[scenario] || !predictionAssets[scenario][y]) {
                throw new Error('Missing predictor asset for ' + scenario + ' ' + y);
            }
        }
    }

    var bands2 = ['predicted_RF', 'predicted_GBT', 'predicted_CART'];
    var bands3 = ['predicted_RF_final', 'predicted_GBT_final', 'predicted_CART_final'];
    var band_P = {min: 0, max: 200,
        palette: ['ffffff','b7f0ae','21f600','0000FF','FDFF92','FF2700','d600ff']};

    for (var s = 0; s < scenarios.length; s++) {
        var scenario = scenarios[s];
        for (var i = startYear; i <= endYear; i += yearStep) {
            var yearOut = String(i);
            var predictors_new = ee.Image(predictionAssets[scenario][i]).select(bands)
                .clip(geometry).toFloat().setDefaultProjection(crs, null, scale);

            // Base models
            var predicted_RF_1 = predictors_new.classify(trained_RF, 'predicted_RF');
            var predicted_GBT_1 = predictors_new.classify(trained_GBT, 'predicted_GBT');
            var predicted_CART_1 = predictors_new.classify(trained_CART, 'predicted_CART');
            var Stacking_predicted_1 = predicted_RF_1.addBands(predicted_GBT_1)
                .addBands(predicted_CART_1).select(bands2);

            // Primary stacking
            var predicted_RF_final_1 = Stacking_predicted_1
                .classify(trained_RF_meta, 'predicted_RF_final');
            var predicted_GBT_final_1 = Stacking_predicted_1
                .classify(trained_GBT_meta, 'predicted_GBT_final');
            var predicted_CART_final_1 = Stacking_predicted_1
                .classify(trained_CART_meta, 'predicted_CART_final');
            var Stacking_predicted_final_1 = predicted_RF_final_1
                .addBands(predicted_GBT_final_1).addBands(predicted_CART_final_1).select(bands3);

            // Secondary stacking: selected RF model
            var predicted_RF_final_final_1 = Stacking_predicted_final_1
                .classify(trained_RF_meta_final, 'carbon_stock')
                .updateMask(vegetationMask).clip(geometry).toFloat().setDefaultProjection(crs, null, scale);
            Map.addLayer(predicted_RF_final_final_1, band_P,
                'Biomass carbon ' + scenario + ' ' + yearOut, false);

            if (exportResults) {
                Export.image.toDrive({
                    image: predicted_RF_final_final_1,
                    folder: exportFolder,
                    description: 'Biomass_carbon_' + scenario + '_' + yearOut,
                    fileNamePrefix: 'Biomass_carbon_' + scenario + '_' + yearOut,
                    region: geometry,
                    crs: crs,
                    scale: scale,
                    maxPixels: 1e13
                });
            }
        }
    }
}

predictBiomassSSP();
