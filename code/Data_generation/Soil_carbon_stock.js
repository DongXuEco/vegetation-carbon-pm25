/*
 * Train SOC or SIC models on 2001-2022 mean predictors, then predict each year.
 * Fit each carbon pool separately; keep predictor bands and units consistent.
 * Output units and soil depth follow the sample target.
 */
// Common valid vegetation domain for all variables, years and scenarios.
// Single-band raster: 1 = included, 0 = excluded; use the same asset in all scripts.
var vegetationMaskAsset = '';

if (!vegetationMaskAsset) {
    throw new Error('Fill in the common vegetation mask asset ID.');
}


var predictorsAsset = ''; // Multiband asset of prepared 2001-2022 mean predictors.
// Point samples contain SOC or SIC stock density for a consistent soil depth.
// Fit one carbon pool per run. Set targetProperty to its numeric field name.
var samplesAsset = '';
var soilCarbonType = 'SOC'; // 'SOC' or 'SIC'.
var targetProperty = soilCarbonType.toLowerCase(); // Change if the field is named differently.
var selectedBands = [];   // Empty: use all predictor bands; exclude target/ID bands.
var predictionAssets = {  // Annual images with the same band names and units.
    2001: '',
    2022: ''
};
var trainingStartYear = 2001;
var trainingEndYear = 2022;
var trainingPeriod = trainingStartYear + '_' + trainingEndYear;
var scale = 10000;
var crs = 'EPSG:4326';
var split = 0.8;
var seed = 2;
var outerFolds = 10;
var innerFolds = 10;       // Cross-fitting within each training subset.
var finalModel = soilCarbonType === 'SOC' ? 'RF' : 'GBT'; // Selected models in the manuscript.
var runValidation = true;
var runAnnualPrediction = false;
var exportResults = true;
var predictionStartYear = 2001;
var predictionEndYear = 2022;
var exportFolder = 'Soil_carbon_' + soilCarbonType;

if (!predictorsAsset || !samplesAsset) {
    throw new Error('Fill in predictorsAsset and samplesAsset.');
}
if (['SOC', 'SIC'].indexOf(soilCarbonType) < 0) {
    throw new Error('soilCarbonType must be SOC or SIC.');
}
if (['RF', 'GBT', 'CART'].indexOf(finalModel) < 0) {
    throw new Error('finalModel must be RF, GBT or CART.');
}
if (runAnnualPrediction) {
    for (var year = predictionStartYear; year <= predictionEndYear; year++) {
        if (!predictionAssets[year]) {
            throw new Error('Missing predictor asset for ' + year);
        }
    }
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

// Ten-fold validation, confined to the development sample set.
var foldMetrics = ee.FeatureCollection([]);
var foldPredictions = ee.FeatureCollection([]);
var foldMaps = [];
if (runValidation) {
    var outerData = assignFolds(trainingData, outerFolds, seed + 10, '_outer_fold');
    print('Development samples by outer fold', outerData.aggregate_histogram('_outer_fold'));
    for (var fold = 0; fold < outerFolds; fold++) {
        var foldTrain = outerData.filter(ee.Filter.neq('_outer_fold', fold));
        var foldTest = outerData.filter(ee.Filter.eq('_outer_fold', fold));
        var fittedFold = fitStack(foldTrain, seed + 100 + fold * 10000);
        var predictedFold = predictStackTable(foldTest, fittedFold);
        foldPredictions = foldPredictions.merge(predictedFold.select(
            ['_sample_id', '_pixel_id', '_outer_fold', targetProperty].concat(outputNames)));
        for (var m = 0; m < outputNames.length; m++) {
            foldMetrics = foldMetrics.merge(ee.FeatureCollection([
                metrics(predictedFold, outputNames[m], 'development_cv', fold)
            ]));
        }
        foldMaps.push(predictStackImage(predictors, fittedFold, finalModel));
    }
    for (var k = 0; k < outputNames.length; k++) {
        foldMetrics = foldMetrics.merge(ee.FeatureCollection([
            metrics(foldPredictions, outputNames[k], 'development_cv_pooled', -1)
        ]));
    }
    var predictionSD = ee.ImageCollection.fromImages(foldMaps)
        .reduce(ee.Reducer.sampleStdDev()).rename('prediction_sd');
    var uncertainty = predictionSD.multiply(1.96).rename('uncertainty_1_96_sd');
    exportImage(predictionSD.addBands(uncertainty), soilCarbonType + '_uncertainty_mean_environment_' + trainingPeriod);
    if (exportResults) {
        Export.table.toDrive({collection: foldMetrics, folder: exportFolder,
            description: soilCarbonType + '_cv_metrics', fileFormat: 'CSV'});
        Export.table.toDrive({collection: foldPredictions, folder: exportFolder,
            description: soilCarbonType + '_cv_predictions', fileFormat: 'CSV',
            selectors: ['_sample_id', '_pixel_id', '_outer_fold', targetProperty].concat(outputNames)});
    }
}

// The independent test set is used only after fitting and model choice.
var fittedFinal = fitStack(trainingData, seed + 1000000);
var testPredictions = predictStackTable(validationData, fittedFinal);
var testMetrics = ee.FeatureCollection([
    metrics(testPredictions, 'secondary_' + finalModel, 'independent_test', -1)
]);
if (exportResults) {
    Export.table.toDrive({collection: testMetrics, folder: exportFolder,
        description: soilCarbonType + '_test_metrics', fileFormat: 'CSV'});
    Export.table.toDrive({collection: testPredictions, folder: exportFolder,
        description: soilCarbonType + '_test_predictions', fileFormat: 'CSV',
        selectors: ['_sample_id', '_pixel_id', targetProperty, 'secondary_' + finalModel]});
}
var predicted_RF_final_final = predictStackImage(predictors, fittedFinal, finalModel);
var band_P = {min: 0, max: 200,
    palette: ['ffffff','b7f0ae','21f600','0000FF','FDFF92','FF2700','d600ff']};
Map.addLayer(predicted_RF_final_final, band_P, soilCarbonType + ' under mean environment ' + trainingPeriod, false);
exportImage(predicted_RF_final_final, soilCarbonType + '_mean_environment_' + trainingPeriod);

// Apply the mean-environment model to every annual predictor image, including 2010.
if (runAnnualPrediction) {
    for (var i = predictionStartYear; i <= predictionEndYear; i++) {
        var predictors_new = ee.Image(predictionAssets[i]).select(bands)
            .clip(geometry).toFloat().setDefaultProjection(crs, null, scale);
        var predicted_RF_final_final_1 = predictStackImage(predictors_new, fittedFinal, finalModel);
        exportImage(predicted_RF_final_final_1, soilCarbonType + '_stock_' + i);
    }
}
