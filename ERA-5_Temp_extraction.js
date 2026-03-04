// Import the FeatureCollection of MVCA polygons
var district = ee.FeatureCollection('projects/era-5-amt/assets/district');

// Print the FeatureCollection to check the properties of a sample feature
var sample_feature = district.first();
print('Sample Feature Properties:', sample_feature.propertyNames());

// Visualize the district polygons on the map
Map.addLayer(district, {color: 'red'}, 'district Polygons');
Map.centerObject(district, 8);

// Define the date range
var start_date = ee.Date('2009-01-01');
var end_date = ee.Date('2010-01-01');

// Load the ERA5 Hourly dataset for temperature
var era5_temp = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
                .select('temperature_2m')
                .filterDate(start_date, end_date);

// Function to extract average temperature for each polygon
function extractERA5Data(feature) {
  // Map over the hourly temperature collection to create a time series
  var temp_ts = era5_temp.map(function(image) {
    var date = image.date().format('YYYY-MM-dd HH:mm');
    
    // Reduce the image over the polygon to get the mean temperature
    var temp_mean = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: feature.geometry(),
      scale: 5000,  // Adjust the scale to match ERA5 resolution (~9km)
      crs: 'EPSG:4326',  // Use the appropriate projection for ERA5
      bestEffort: true,  // Use this to avoid computation limits
      maxPixels: 1e9  // Increase maxPixels if needed
    }).get('temperature_2m');  // Get the mean temperature for the polygon
    
     // Check if temp_mean is null and handle accordingly
    var temp_celsius = ee.Algorithms.If(
      ee.Algorithms.IsEqual(temp_mean, null),
      null,  // Return null if temp_mean is null
      ee.Number(temp_mean).subtract(273.15)  // Convert to Celsius
    );
    
    // Get the correct property names from the feature
    var site_id = feature.get('Name');  // Replace 'site' with the correct property name if needed
    //var id_num = feature.get('id_num');  // Get the 'id_num' property
    
    // Return a new feature with the extracted temperature value
    return ee.Feature(null, {
      'date': date,
      'temperature': temp_celsius,
      'site_id': site_id,
     // 'id_num': id_num,  // Include the 'id_num' property
    //  'latitude': feature.geometry().centroid().coordinates().get(1),  // Centroid for reference
    //  'longitude': feature.geometry().centroid().coordinates().get(0)  // Centroid for reference
    });
  });

  return temp_ts;
}

// Apply the function to each polygon in the collection
var results = district.map(extractERA5Data).flatten();

// Print results to the console
print('Results:', results);

// Export the results to a CSV file
Export.table.toDrive({
  collection: results,
  description: 'ERA5_Temperature_Extraction_Hourly_district',
  fileFormat: 'CSV',
  selectors: ['date', 'temperature', 'site_id']  
});
