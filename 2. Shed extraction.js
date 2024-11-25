// ================== 1. 加载研究区域和样本数据 ==================

// 加载 FAO GAUL 2015 一级行政区划数据
var gaul = ee.FeatureCollection('FAO/GAUL/2015/level1');

// 过滤出中国的行政区划
var china_adm1 = gaul.filter(ee.Filter.eq('ADM0_NAME', 'China'));

// 定义感兴趣的省份列表
var provinces_of_interest = [
  'Hunan Sheng', 
  'Jiangxi Sheng', 
  'Guangdong Sheng', 
  'Fujian Sheng', 
  'Guangxi Zhuangzu Zizhiqu', 
  'Guizhou Sheng'
];

// 过滤出感兴趣的省份并合并为一个几何对象
var roi = china_adm1.filter(ee.Filter.inList('ADM1_NAME', provinces_of_interest))
  .geometry();

Map.centerObject(roi, 6);
Map.addLayer(roi, {color: 'blue'}, '研究区域');

// 加载样本数据
var greenhouseSamples2003 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/GreenhouseSamples_2003");
var greenhouseSamples2013 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/GreenhouseSamples_2013");
var greenhouseSamples2023 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/GreenhouseSamples_2023");
var nonGreenhouseSamples2003 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/NonGreenhouseSamples_2003");
var nonGreenhouseSamples2013 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/NonGreenhouseSamples_2013");
var nonGreenhouseSamples2023 = ee.FeatureCollection("projects/ee-1261423515/assets/DP/NonGreenhouseSamples_2023");

// ================== 2. 定义 Landsat 数据处理函数 ==================

// 定义云掩膜和预处理函数
function maskLandsat(image) {
  var qa = image.select('QA_PIXEL');
  var cloudShadowBitMask = (1 << 4);
  var snowBitMask = (1 << 5);
  var cloudBitMask = (1 << 3);
  var cirrusBitMask = (1 << 2);

  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(snowBitMask).eq(0))
               .and(qa.bitwiseAnd(cloudBitMask).eq(0))
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  
  return image.updateMask(mask).select('SR_B.*')
              .multiply(0.0000275).add(-0.2)
              .copyProperties(image, ['system:time_start']);
}

// 加载并预处理 Landsat 影像的函数
function loadLandsatImage(year, roi) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);

  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(maskLandsat)
    .select(
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
    );

  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(maskLandsat)
    .select(
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
    );

  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .map(maskLandsat)
    .select(
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
      ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
    );

  // 根据年份选择合适的影像集合
  var imageCollection;
  if (year <= 2011) {
    imageCollection = l5.merge(l7);
  } else if (year <= 2013) {
    imageCollection = l5.merge(l7);
  } else if (year <= 2019) {
    imageCollection = l7.merge(l8);
  } else {
    imageCollection = l8;
  }

  // 生成年度中值合成影像
  return imageCollection.median().clip(roi);
}

// 添加指数的函数
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI');
  var ndwi = image.normalizedDifference(['B3', 'B5']).rename('NDWI');
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    }).rename('EVI');

  var indices = [ndvi, ndwi, evi];

  return image.addBands(indices);
}

// ================== 3. 计算纹理特征的函数 ==================

function addTextureFeatures(image) {
  // 计算灰度图像
  var gray = image.expression(
    '(0.3 * NIR) + (0.59 * RED) + (0.11 * GREEN)', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'GREEN': image.select('B3')
  }).rename('gray');

  // 计算 GLCM 纹理特征
  var glcm = gray.multiply(100).toInt().glcmTexture({size: 3});

  var contrast = glcm.select('gray_contrast');  // 对比度
  var dissimilarity = glcm.select('gray_diss'); // 差异性
  var entropy = glcm.select('gray_ent');    // 熵
  var angularSecondMoment = glcm.select('gray_asm'); // 角二阶矩

  return image.addBands([contrast, dissimilarity, entropy, angularSecondMoment]);
}

// ================== 4. 面向对象的图像分割函数 ==================

function segmentImage(image) {
  // 设置种子
  var seeds = ee.Algorithms.Image.Segmentation.seedGrid(30);

  // 使用 SNIC 算法进行图像分割
  var snic = ee.Algorithms.Image.Segmentation.SNIC({
    image: image,
    size: 32,
    compactness: 3,
    connectivity: 8,
    neighborhoodSize: 128,
    seeds: seeds
  });

  // 获取分割后的聚类标签
  var clusters = snic.select('clusters');

  // 计算每个对象的统计特征（均值）
  var meanImage = image.addBands(clusters).reduceConnectedComponents({
    reducer: ee.Reducer.mean(),
    labelBand: 'clusters'
  });

  // 计算每个对象的标准差
  var stdDevImage = image.addBands(clusters).reduceConnectedComponents({
    reducer: ee.Reducer.stdDev(),
    labelBand: 'clusters'
  });

  // 合并对象特征
  return meanImage.addBands(stdDevImage);
}

// ================== 5. 对每一年进行处理和分类 ==================

var years = [2003, 2013, 2023];

years.forEach(function(year) {
  print('Processing year:', year);

  // 加载对应年份的样本数据
  var greenhouseSamples, nonGreenhouseSamples;
  if (year === 2003) {
    greenhouseSamples = greenhouseSamples2003;
    nonGreenhouseSamples = nonGreenhouseSamples2003;
  } else if (year === 2013) {
    greenhouseSamples = greenhouseSamples2013;
    nonGreenhouseSamples = nonGreenhouseSamples2013;
  } else if (year === 2023) {
    greenhouseSamples = greenhouseSamples2023;
    nonGreenhouseSamples = nonGreenhouseSamples2023;
  }

  // 为样本添加类别标签
  function setLabel(feature, label) {
    return feature.set({'Landcover': label});
  }
  var greenhouseSamplesLabeled = greenhouseSamples.map(function(feature) {
    return setLabel(feature, 1); // 大棚类别为 1
  });
  var nonGreenhouseSamplesLabeled = nonGreenhouseSamples.map(function(feature) {
    return setLabel(feature, 0); // 非大棚类别为 0
  });

  // 合并样本
  var combinedTrainingSet = greenhouseSamplesLabeled.merge(nonGreenhouseSamplesLabeled);

  // 可视化训练样本（可选）
  var visualizationParams1 = { color: 'FF0000', pointRadius: 2 }; // 红色，大棚
  var visualizationParams0 = { color: '0000FF', pointRadius: 2 }; // 蓝色，非大棚
  Map.addLayer(greenhouseSamplesLabeled, visualizationParams1, 'Greenhouse Samples ' + year);
  Map.addLayer(nonGreenhouseSamplesLabeled, visualizationParams0, 'Non-Greenhouse Samples ' + year);

  // 加载并预处理对应年份的 Landsat 影像
  var image = loadLandsatImage(year, roi);
  image = addIndices(image); // 添加指数
  image = addTextureFeatures(image); // 添加纹理特征

  // 对影像进行分割，得到对象特征
  var segmentedImage = segmentImage(image);

  // 合并影像和对象特征
  var data = image.addBands(segmentedImage);

  // 选择用于分类的波段
  var bands = data.bandNames();

  // 添加随机数字段
  var sampleData = combinedTrainingSet.randomColumn('random');

  // 随机拆分样本点为训练样本和验证样本
  var sample_training = sampleData.filter(ee.Filter.lte('random', 0.8));
  var sample_validate  = sampleData.filter(ee.Filter.gt('random', 0.8));

  // 利用样本点拾取特征值用于模型训练和验证
  var training = data.select(bands).sampleRegions({
    collection: sample_training,
    properties: ['Landcover'],
    scale: 30,
    tileScale: 16
  });

  var validation = data.select(bands).sampleRegions({
    collection: sample_validate,
    properties: ['Landcover'],
    scale: 30,
    tileScale: 16
  });

  // 训练随机森林分类器
  var rf = ee.Classifier.smileRandomForest({
    numberOfTrees: 100,
    bagFraction: 0.8
  }).train({
    features: training,
    classProperty: 'Landcover',
    inputProperties: bands
  });

  // 对影像进行分类
  var classified = data.select(bands).classify(rf);

  // 验证分类器
  var validated = validation.classify(rf);
  var testAccuracy = validated.errorMatrix('Landcover', 'classification');
  print('Error matrix for year ' + year + ':', testAccuracy);
  print('Overall accuracy for year ' + year + ':', testAccuracy.accuracy());

  // 导出分类结果
  Export.image.toDrive({
    image: classified.clip(roi).uint8(),
    description: 'Greenhouse_Classification_' + year,
    folder: 'Greenhouse_Classification',
    fileNamePrefix: 'Greenhouse_Classification_' + year,
    region: roi,
    scale: 30,
    maxPixels: 1e13,
    crs: 'EPSG:4326',
    formatOptions: {
      cloudOptimized: true
    }
  });

  // 导出分类器信息（可选）
  Export.table.toDrive({
    collection: ee.FeatureCollection([ee.Feature(null, {
      'year': year,
      'trainingAccuracy': rf.confusionMatrix().accuracy(),
      'validationAccuracy': testAccuracy.accuracy()
    })]),
    description: 'Classification_Accuracy_' + year,
    folder: 'Greenhouse_Classification',
    fileNamePrefix: 'Classification_Accuracy_' + year,
    fileFormat: 'CSV'
  });

});

