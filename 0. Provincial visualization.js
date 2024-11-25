// 加载 FAO GAUL 2015 一级行政区划数据
var gaul = ee.FeatureCollection('FAO/GAUL/2015/level1');

// 过滤出中国的行政区划
var china_adm1 = gaul.filter(ee.Filter.eq('ADM0_NAME', 'China'));

// 打印中国所有省份的名称
var province_names = china_adm1.aggregate_array('ADM1_NAME');
print('中国的省份名称：', province_names);

// 定义感兴趣的省份列表（与数据集名称一致）
var provinces_of_interest = [
  'Hunan Sheng', 
  'Jiangxi Sheng', 
  'Guangdong Sheng', 
  'Fujian Sheng', 
  'Guangxi Zhuangzu Zizhiqu', 
  'Guizhou Sheng'
];

// 过滤出感兴趣的省份
var selected_provinces = china_adm1.filter(ee.Filter.inList('ADM1_NAME', provinces_of_interest));

// 可视化设置
var provinceStyle = {
  color: 'blue',
  fillColor: '00000000',
  width: 2
};

// 将选定的省份添加到地图
Map.centerObject(selected_provinces, 5);
Map.addLayer(selected_provinces.style(provinceStyle), {}, '选定省份');

// 可选：显示所有省份用于参考
Map.addLayer(china_adm1.style({color: 'gray', fillColor: '00000000', width: 1}), {}, '中国所有省份');
