var palettes = require('users/gena/packages:palettes');
var ui2 = require('users/gena/packages:ui');

// Load Landsat 5 and Landsat 8 images

var landsat5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_TOA')
    .filterBounds(table)
    .filterDate('1987-01-01', '1987-03-31')
    .mean();

var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterBounds(table)
    .filterDate('2024-01-01', '2024-03-31')
    .mean();

// Calculate NIR/SWIR ratio
var nir_swir_ratio_l5 = landsat5.select('B4').divide(landsat5.select('B5'));
var nir_swir_ratio_l8 = landsat8.select('B5').divide(landsat8.select('B6'));

// Display the result
//Map.addLayer(nir_swir_ratio_l5, {min: 0, max: 80, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 5');
//Map.addLayer(nir_swir_ratio_l8, {min: 0, max: 80, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 8');

var HSpalette = palettes.cmocean.Ice[7];
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

// Create the satellite basemap
var satelliteBasemap = ee.Image('srtm90_v4');

// Calculate NIR/SWIR ratio
var nir_swir_ratio_l5_layer = ui.Map.Layer(nir_swir_ratio_l5.mask(nir_swir_ratio_l5.gt(3)), {min: 0, max: 80, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 5');
var nir_swir_ratio_l8_layer = ui.Map.Layer(nir_swir_ratio_l8.mask(nir_swir_ratio_l8.gt(3)), {min: 0, max:80, palette: ['black', 'cyan', 'yellow', 'white']}, 'NIR/SWIR Ratio Landsat 8');

var leftMap = ui.Map();
leftMap.add(satelliteBasemap, {}, 'Satellite Basemap'); // Add satellite basemap
leftMap.addLayer(nir_swir_ratio_l5.mask(nir_swir_ratio_l5.gt(3)), {min: 0, max: 80, palette: HSpalette}, 'NIR/SWIR Ratio Landsat 5');
var palettes = require('users/gena/packages:palettes');

// Create the right map panel.
var rightMap = ui.Map();
rightMap.add(satelliteBasemap, {}, 'Satellite Basemap'); // Add satellite basemap
rightMap.addLayer(nir_swir_ratio_l8.mask(nir_swir_ratio_l8.gt(3)), {min: 0, max: 80, palette: HSpalette}, 'NIR/SWIR Ratio Landsat 8');
rightMap.setControlVisibility({zoomControl: false});

// Create a label for the right map panel.
var rightLabel = ui.Label('Detección de glaciares con Landsat 8 2024');
rightMap.add(rightLabel);

var rileftLabel = ui.Label('Detección de glaciares con Landsat 5 1987');
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


var vis = {min: 0, max: 80, palette: 'ice'};


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

// Add the legendPanel to the map.
var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
leftMap.add(legendPanel);

add original landsat 5 and 8 images to the map behind every other Layer
leftMap.addLayer(landsat5, {bands: ['B3', 'B2', 'B1'], min: 0, max: 0.3}, 'Landsat 5');

