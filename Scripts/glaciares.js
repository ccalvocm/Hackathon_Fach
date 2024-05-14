// Load Landsat 5 and Landsat 8 images

var landsat5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_TOA')
    .filterBounds(table)
    .filterDate('1987-01-01', '1987-03-31')
    .mean();

var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterBounds(table)
    .filterDate('2024-01-01', '2024-03-31')
    .mean();
    
var landsat82 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterBounds(table)
    .filterDate('2022-01-01', '2022-02-1')
    .mean();    

// Calculate NIR/SWIR ratio
var nir_swir_ratio_l5 = landsat5.select('B4').divide(landsat5.select('B5'));
var nir_swir_ratio_l8 = landsat8.select('B5').divide(landsat8.select('B6'));

// export to Google Drive nir_swir_ratio_l5 and nir_swir_ratio_l8
Export.image.toDrive({
    image: nir_swir_ratio_l5,
    description: 'nir_swir_ratio_l5',
    scale: 30,
    region: table.geometry().bounds(),
    maxPixels: 1e13,
    crs: 'EPSG:32719'
});

Export.image.toDrive({
    image: nir_swir_ratio_l8,
    description: 'nir_swir_ratio_l8',
    scale: 30,
    region: table.geometry().bounds(),
    maxPixels: 1e13,
    crs: 'EPSG:32719'
});


var palettes = require('users/gena/packages:palettes');
var HSpalette = palettes.kovesi.linear_blue_5_95_c73[7];
function makeColorBarParams(palette) {
  return {
    bbox: [0, 0, 7, 1],
    dimensions: '100x10',
    format: 'png',
    min: 0,
    max: 7,
    palette: HSpalette,
  };
}
// Calculate NIR/SWIR ratio
var nir_swir_ratio_l5_layer = ui.Map.Layer(nir_swir_ratio_l5.mask(nir_swir_ratio_l5.gt(4)), {min: 0, max: 20, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 5');
var nir_swir_ratio_l8_layer = ui.Map.Layer(nir_swir_ratio_l8.mask(nir_swir_ratio_l8.gt(4)), {min: 0, max:20, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 8');

var leftMap = ui.Map();
leftMap.addLayer(landsat5, {bands: ['B3', 'B2', 'B1'], min: 0, max: 0.3}, 'Landsat 8');
leftMap.addLayer(nir_swir_ratio_l5.mask(nir_swir_ratio_l5.gt(4)), {min: 0, max: 20, palette: HSpalette}, 'NIR/SWIR Ratio Landsat 5');
leftMap.setControlVisibility({zoomControl: true});

// Create the right map panel.
var rightMap = ui.Map();
rightMap.addLayer(landsat82, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3}, 'Landsat 8');
rightMap.addLayer(nir_swir_ratio_l8.mask(nir_swir_ratio_l8.gt(4)), {min: 0, max: 20, palette: HSpalette}, 'NIR/SWIR Ratio Landsat 8');
rightMap.setControlVisibility({zoomControl: false});


// Create a label for the right map panel.
var rightLabel = ui.Label('Detección de glaciares con Landsat 8 2024', {fontWeight: 'bold'});
rightMap.add(rightLabel);

var rileftLabel = ui.Label('Detección de glaciares con Landsat 5 1987', {fontWeight: 'bold'});
leftMap.add(rileftLabel);

// Create a SplitPanel to hold the two maps.
var splitPanel = ui.SplitPanel({
  firstPanel: leftMap,
  secondPanel: rightMap,
  orientation: 'horizontal',
  wipe: true
});

// Set the SplitPanel as the only thing in the UI root.
ui.root.widgets().reset([splitPanel]);

// Link the two maps.
var linker = ui.Map.Linker([leftMap, rightMap]);

var ui2 = require('users/gena/packages:ui');

var vis = {min: 3, max: 20, palette: 'ice'};


// Create the color bar for the legend.
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0),
  params: makeColorBarParams(vis.palette),
  style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
});

var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0),
  params: makeColorBarParams(vis.palette),
  style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
});

// Create a panel with three numbers for the legend.
var legendLabels = ui.Panel({
  widgets: [
    ui.Label(vis.min, {margin: '4px 8px'}),
    ui.Label(
        ((vis.max-vis.min) / 2+vis.min),
        {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(vis.max, {margin: '4px 8px'})
  ],
  layout: ui.Panel.Layout.flow('horizontal')
});

var legendTitle = ui.Label({
    value: 'Índices NIR/SWIR para la detección de glaciares',
    style: {fontWeight: 'bold'}
});

// Add the legendPanel to the map displaced to the left
var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
leftMap.add(legendPanel)
        .setControlVisibility({position: 'bottomleft'})
        .setControlVisibility({position: 'topleft'}); // Add this line to set the legendPanel to the left side of the map.
