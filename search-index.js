var searchIndex = JSON.parse('{\
"solar":{"doc":"","t":[6,17,3,3,17,4,3,13,13,13,4,13,12,5,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,11,11,0,0,12,11,12,11,11,11,11,11,11,11,11,11,11,11,12,12,3,4,13,13,17,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,3,17,12,11,11,5,11,11,11,5,11,11,11,12,11,12,12,11,5,11,11,5,12,11,11,11,11,11],"n":["Float","PI","PerezSky","ReinhartSky","SOLAR_CONSTANT","SkyUnits","Solar","Solar","Solar","Standard","Time","Visible","acc_bins","air_mass","b","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","declination","equation_of_time","from","from","hour_angle","into","into","latitude","longitude","mf","n_bins","n_rows","new","normal_extraterrestrial_radiation","perez","reinhart_sky","row_max_sin","solar_standard_time_difference","standard_meridian","sun_position","sunrise_sunset","to_owned","try_from","try_from","try_into","try_into","type_id","type_id","unwrap_solar_time","unwrap_standard_time","0","0","PerezSky","SkyUnits","Solar","Visible","WHTEFFICACY","borrow","borrow","borrow_mut","borrow_mut","calc_params","clearness_category","clone","clone_into","diffuse_illuminance_ratio","direct_illuminance_ratio","from","from","gen_sky_vec","get_sky_func_standard_time","into","into","precipitable_water_content","sky_brightness","sky_clearness","to_owned","try_from","try_from","try_into","try_into","type_id","type_id","update_sky_vec","ReinhartSky","TNAZ","acc_bins","bin_dir","bin_solid_angle","bins_in_row","bins_row","borrow","borrow_mut","cone_solid_angle","dir_to_bin","from","into","mf","n_bins","n_bins","n_rows","new","patch_solid_angle","raccum","row_altitude","row_height","row_max_sin","sin_altitude_to_row","try_from","try_into","type_id","xy_to_bin"],"q":["solar","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","solar::Time","","solar::perez","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","solar::reinhart_sky","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"d":["","","This struct packs functions from two sources:","","W/m2","Specifies which units do we want returned from the sky …","Solar calculations library. Based on Duffie and Beckman’…","Solar Radiance (i.e., W/m2.sr)","Time is in Solar time","Time is in Standard time","Solar or Standard time, containing the day of the year ‘n…","Visible Radiance (i.e., W/m2.sr, but only visible)","An array with the number of bins accumulated before each …","Calculates the Air-mass .PerezSky","Equation 1.4.2 in the Book.","","","","","","","Declination (in Radians), according to Equation 1.6.1B","The Equation of Time based on the day of year (can have …","","","Returns the hour angle in degrees    ","","","Latitude in Radians","Longitude (in Radians)","Subdivition scheme","The number of elements in this sky discretization …","The number of rows, including Cap","Builds a Solar object from  a Latitude, Longitude and …","Normal extraterrestrial radiation (Gon) Equation 1.4.1b …","","","An array with the Sin(max angle) for every row.","Returns the difference between the solar and the standard …","Standard meridian (in Radians)","Builds a vector that points towards the sun.","Gets the sunset time (equation 1.6.10) n should be in …","","","","","","","","Returns the content of a Time enum. Transforms to Solar if …","Returns the content of a Time enum. Transforms to Standard …","","","This struct packs functions from two sources:","Specifies which units do we want returned from the sky …","Solar Radiance (i.e., W/m2.sr)","Visible Radiance (i.e., W/m2.sr, but only visible)","","","","","","Calculates the Perez sky params","Get the sky clearness category (Table 1 in the 1990 paper)","","","Calculates the diffuse illuminance/diffuse radiance ratio …","","","","","Returns a function that can be used to calculate the Perez…","","","Equation 3… dew_point_temp in C","Equation 2","Equation 1","","","","","","","","","","The number of patches in teach of Tregenza’s (i.e., …","An array with the number of bins accumulated before each …","Gets a <code>Vector3D</code> pointing to the centre the bin.","Calculates te solid angle of a certain bin","Returns the number of elements in row <code>row</code>, starting from 0.","Calculates the row in which a certain bin is located","","","Calculates the solid angle of a Cone.","Returns the bin number pointed by direction <code>dir</code>","","","Subdivition scheme","Calculates the number of total bins in a Reinhart’s …","The number of elements in this sky discretization …","The number of rows, including Cap","Creates  a new Reinhart sky discretization","","The number of sky elements accumulated up to row <code>row</code> …","Returns the minimum and maximum altitude of a row, in …","Returns the height of a full row (the cap is only half …","An array with the Sin(max angle) for every row.","Gets the sin of an altitude","","","","Gets the position of the bin within a row."],"i":[0,0,0,0,0,0,0,1,2,2,0,1,3,0,4,4,2,4,2,2,2,4,4,4,2,4,4,2,4,4,3,3,3,4,4,0,0,3,4,4,4,4,2,4,2,4,2,4,2,4,4,5,6,0,0,1,1,0,7,1,7,1,7,7,1,1,7,7,7,1,7,7,7,1,7,7,7,1,7,1,7,1,7,1,7,0,0,3,3,3,0,3,3,3,0,3,3,3,3,3,3,3,3,0,3,3,0,3,3,3,3,3,3],"f":[null,null,null,null,null,null,null,null,null,null,null,null,null,[[["f64",15]],["f64",15]],[[["f64",15]],["f64",15]],[[]],[[]],[[]],[[]],[[],["time",4]],[[]],[[["f64",15]],["f64",15]],[[["f64",15]],["f64",15]],[[]],[[]],[[["time",4]],["f64",15]],[[]],[[]],null,null,null,null,null,[[["f64",15],["f64",15],["f64",15]]],[[["f64",15]],["f64",15]],null,null,null,[[["f64",15]],["f64",15]],null,[[["time",4]],["option",4,[["vector3d",3]]]],[[["f64",15]]],[[]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],[[["time",4]],["f64",15]],[[["time",4]],["f64",15]],null,null,null,null,null,null,null,[[]],[[]],[[]],[[]],[[["f64",15],["f64",15],["f64",15]]],[[["f64",15]],["usize",15]],[[],["skyunits",4]],[[]],[[["f64",15],["f64",15],["f64",15],["usize",15]],["f64",15]],[[["f64",15],["f64",15],["f64",15],["usize",15]],["f64",15]],[[]],[[]],[[["usize",15],["solar",3],["date",3],["currentweather",3],["skyunits",4],["f64",15],["bool",15],["bool",15]],["result",4,[["matrix",6],["string",3]]]],[[["skyunits",4],["solar",3],["date",3],["f64",15],["f64",15],["f64",15]],["box",3,[["fn",8]]]],[[]],[[]],[[["f64",15]],["f64",15]],[[["f64",15],["f64",15],["f64",15]],["f64",15]],[[["f64",15],["f64",15],["f64",15]],["f64",15]],[[]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],[[["matrix",6],["usize",15],["solar",3],["date",3],["currentweather",3],["skyunits",4],["f64",15],["bool",15],["bool",15]],["result",4,[["string",3]]]],null,null,null,[[["usize",15]],["vector3d",3]],[[["usize",15]],["f64",15]],[[["usize",15],["usize",15]],["usize",15]],[[["usize",15]],["usize",15]],[[]],[[]],[[["f64",15]],["f64",15]],[[["vector3d",3]],["usize",15]],[[]],[[]],null,[[["usize",15]],["usize",15]],null,null,[[["usize",15]]],[[["f64",15],["f64",15],["usize",15]],["f64",15]],[[["usize",15]],["usize",15]],[[["usize",15]]],[[["usize",15]],["f64",15]],null,[[["f64",15]],["usize",15]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[["usize",15],["f64",15],["f64",15]],["usize",15]]],"p":[[4,"SkyUnits"],[4,"Time"],[3,"ReinhartSky"],[3,"Solar"],[13,"Solar"],[13,"Standard"],[3,"PerezSky"]]}\
}');
if (window.initSearch) {window.initSearch(searchIndex)};