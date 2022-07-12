var searchIndex = JSON.parse('{\
"solar":{"doc":"","t":[6,17,3,3,17,4,3,13,13,13,4,13,12,5,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,11,11,0,0,12,11,12,11,11,11,11,11,11,11,11,11,11,11,12,12,3,4,13,13,17,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,3,17,12,11,11,5,11,11,11,5,11,11,11,12,11,12,12,11,5,11,11,5,12,11,11,11,11,11],"n":["Float","PI","PerezSky","ReinhartSky","SOLAR_CONSTANT","SkyUnits","Solar","Solar","Solar","Standard","Time","Visible","acc_bins","air_mass","b","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","declination","equation_of_time","from","from","hour_angle","into","into","latitude","longitude","mf","n_bins","n_rows","new","normal_extraterrestrial_radiation","perez","reinhart_sky","row_max_sin","solar_standard_time_difference","standard_meridian","sun_position","sunrise_sunset","to_owned","try_from","try_from","try_into","try_into","type_id","type_id","unwrap_solar_time","unwrap_standard_time","0","0","PerezSky","SkyUnits","Solar","Visible","WHTEFFICACY","borrow","borrow","borrow_mut","borrow_mut","calc_params","clearness_category","clone","clone_into","diffuse_illuminance_ratio","direct_illuminance_ratio","from","from","gen_sky_vec","get_sky_func_standard_time","into","into","precipitable_water_content","sky_brightness","sky_clearness","to_owned","try_from","try_from","try_into","try_into","type_id","type_id","update_sky_vec","ReinhartSky","TNAZ","acc_bins","bin_dir","bin_solid_angle","bins_in_row","bins_row","borrow","borrow_mut","cone_solid_angle","dir_to_bin","from","into","mf","n_bins","n_bins","n_rows","new","patch_solid_angle","raccum","row_altitude","row_height","row_max_sin","sin_altitude_to_row","try_from","try_into","type_id","xy_to_bin"],"q":["solar","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","solar::Time","","solar::perez","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","solar::reinhart_sky","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"d":["","","This struct packs functions from two sources:","","W/m2","Specifies which units do we want returned from the sky …","Solar calculations library. Based on Duffie and Beckman’…","Solar Radiance (i.e., W/m2.sr)","Time is in Solar time","Time is in Standard time","Solar or Standard time, containing the day of the year ‘n…","Visible Radiance (i.e., W/m2.sr, but only visible)","An array with the number of bins accumulated before each …","Calculates the Air-mass .PerezSky","Equation 1.4.2 in the Book.","","","","","","","Declination (in Radians), according to Equation 1.6.1B","The Equation of Time based on the day of year (can have …","Returns the argument unchanged.","Returns the argument unchanged.","Returns the hour angle in degrees    ","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Latitude in Radians","Longitude (in Radians)","Subdivition scheme","The number of elements in this sky discretization …","The number of rows, including Cap","Builds a Solar object from  a Latitude, Longitude and …","Normal extraterrestrial radiation (Gon) Equation 1.4.1b …","","","An array with the Sin(max angle) for every row.","Returns the difference between the solar and the standard …","Standard meridian (in Radians)","Builds a vector that points towards the sun.","Gets the sunset time (equation 1.6.10) n should be in …","","","","","","","","Returns the content of a Time enum. Transforms to Solar if …","Returns the content of a Time enum. Transforms to Standard …","","","This struct packs functions from two sources:","Specifies which units do we want returned from the sky …","Solar Radiance (i.e., W/m2.sr)","Visible Radiance (i.e., W/m2.sr, but only visible)","","","","","","Calculates the Perez sky params","Get the sky clearness category (Table 1 in the 1990 paper)","","","Calculates the diffuse illuminance/diffuse radiance ratio …","","Returns the argument unchanged.","Returns the argument unchanged.","","Returns a function that can be used to calculate the Perez…","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Equation 3… dew_point_temp in C","Equation 2","Equation 1","","","","","","","","","","The number of patches in teach of Tregenza’s (i.e., …","An array with the number of bins accumulated before each …","Gets a <code>Vector3D</code> pointing to the centre the bin.","Calculates te solid angle of a certain bin","Returns the number of elements in row <code>row</code>, starting from 0.","Calculates the row in which a certain bin is located","","","Calculates the solid angle of a Cone.","Returns the bin number pointed by direction <code>dir</code>","Returns the argument unchanged.","Calls <code>U::from(self)</code>.","Subdivition scheme","Calculates the number of total bins in a Reinhart’s …","The number of elements in this sky discretization …","The number of rows, including Cap","Creates  a new Reinhart sky discretization","","The number of sky elements accumulated up to row <code>row</code> …","Returns the minimum and maximum altitude of a row, in …","Returns the height of a full row (the cap is only half …","An array with the Sin(max angle) for every row.","Gets the sin of an altitude","","","","Gets the position of the bin within a row."],"i":[0,0,0,0,0,0,0,9,3,3,0,9,17,0,2,2,3,2,3,3,3,2,2,2,3,2,2,3,2,2,17,17,17,2,2,0,0,17,2,2,2,2,3,2,3,2,3,2,3,2,2,18,19,0,0,9,9,0,20,9,20,9,20,20,9,9,20,20,20,9,20,20,20,9,20,20,20,9,20,9,20,9,20,9,20,0,0,17,17,17,0,17,17,17,0,17,17,17,17,17,17,17,17,0,17,17,0,17,17,17,17,17,17],"f":[0,0,0,0,0,0,0,0,0,0,0,0,0,[1,1],[[2,1],1],[[]],[[]],[[]],[[]],[3,3],[[]],[[2,1],1],[[2,1],1],[[]],[[]],[[2,3],1],[[]],[[]],0,0,0,0,0,[[1,1,1],2],[[2,1],1],0,0,0,[[2,1],1],0,[[2,3],[[5,[4]]]],[[2,1]],[[]],[[],6],[[],6],[[],6],[[],6],[[],7],[[],7],[[2,3],1],[[2,3],1],0,0,0,0,0,0,0,[[]],[[]],[[]],[[]],[[1,1,1]],[1,8],[9,9],[[]],[[1,1,1,8],1],[[1,1,1,8],1],[[]],[[]],[[8,2,10,11,9,1,12,12],[[6,[13,14]]]],[[9,2,10,1,1,1],[[16,[15]]]],[[]],[[]],[1,1],[[1,1,1],1],[[1,1,1],1],[[]],[[],6],[[],6],[[],6],[[],6],[[],7],[[],7],[[13,8,2,10,11,9,1,12,12],[[6,[14]]]],0,0,0,[[17,8],4],[[17,8],1],[[8,8],8],[[17,8],8],[[]],[[]],[1,1],[[17,4],8],[[]],[[]],0,[8,8],0,0,[8,17],[[1,1,8],1],[[17,8],8],[[17,8]],[8,1],0,[[17,1],8],[[],6],[[],6],[[],7],[[17,8,1,1],8]],"p":[[15,"f64"],[3,"Solar"],[4,"Time"],[3,"Vector3D"],[4,"Option"],[4,"Result"],[3,"TypeId"],[15,"usize"],[4,"SkyUnits"],[3,"Date"],[3,"CurrentWeather"],[15,"bool"],[6,"Matrix"],[3,"String"],[8,"Fn"],[3,"Box"],[3,"ReinhartSky"],[13,"Solar"],[13,"Standard"],[3,"PerezSky"]]}\
}');
if (typeof window !== 'undefined' && window.initSearch) {window.initSearch(searchIndex)};
if (typeof exports !== 'undefined') {exports.searchIndex = searchIndex};
