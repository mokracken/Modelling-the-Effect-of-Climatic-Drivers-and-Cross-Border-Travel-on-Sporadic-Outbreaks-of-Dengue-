// Import the FeatureCollection of district polygons
var district = ee.FeatureCollection('projects/era-5-amt/assets/district');

// Print the FeatureCollection to check the properties of a sample feature
var sample_feature = district.first();
print('Sample Feature Properties:', sample_feature.propertyNames());

// Visualize the district polygons on the map
Map.addLayer(district, {color: 'blue'}, 'District Polygons');
Map.centerObject(district, 8);

// Define the date range for daily precipitation
var start_date = ee.Date('2009-01-01');
var end_date = ee.Date('2025-01-01');  // For a full year

// Load the CHIRPS Daily dataset for precipitation
var chirps_precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
                .select('precipitation')
                .filterDate(start_date, end_date);

// Function to extract daily precipitation for each polygon
function extractCHIRPSData(feature) {
  var site_id = feature.get('Name');  // District Name

  // Map over the daily precipitation collection to create a time series
  var precip_ts = chirps_precip.map(function(image) {
    var date = image.date().format('YYYY-MM-dd');

    // Reduce the image over the polygon to get the average precipitation
    var precip_avg = image.reduceRegion({
      reducer: ee.Reducer.mean(),  // Use mean instead of sum to avoid extreme values
      geometry: feature.geometry(),
      scale: 5500,  // CHIRPS resolution is ~5.6 km
      crs: 'EPSG:4326',
      bestEffort: true,
      maxPixels: 1e9
    }).get('precipitation');  // Get the average precipitation for the district

    return ee.Feature(null, {
      'date': date,
      'precipitation_mm': precip_avg,  // Already in mm
      'site_id': site_id
    });
  });

  return precip_ts;
}

// Apply the function to each polygon in the collection
var results = district.map(extractCHIRPSData).flatten();

// Print results to the console
print('CHIRPS Precipitation Data:', results);

// Export the results to a CSV file
Export.table.toDrive({
  collection: results,
  description: 'CHIRPS_Precipitation_Mauritius_mm',
  fileFormat: 'CSV',
  selectors: ['date', 'precipitation_mm', 'site_id']
});
