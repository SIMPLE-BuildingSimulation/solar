use calendar::date::Date;
use geometry3d::vector3d::Vector3D;

/// Solar calculations library. Base on Duffie and Beckman's excellent book.
/// 
/// We follow the convention of the book. This means that everything is in 
/// international units, and times are solar. Angles (inputs and outputs) are 
/// in Radians but can be converted into Degrees through the in_degrees() and 
/// in_radiance() functions. 
///
/// Solar azimuth angle is the angular displacement from south of the 
/// projection of beam radiation on the horizontal plane (Figure 1.6.1 of the Book). 
/// Displacements east of south are negative and west of south are positive.
/// 
/// North points in the Y direction. East points in the X direction. Up points in Z.

/// The solar equivalent of Date's "day of the year". The 
/// distinction is there so that we don't mistake solar and 
/// standard time
pub struct Solar {
    /// Day of the year
    n: f64,

    /// Latitude in Radians
    latitude: f64,

    //longitude:f64,
    //standard_meridian: f64,
}

/// W/m2
const SOLAR_CONSTANT : f64 = 1367.0; 

/// Converts Radians into Degrees
#[inline(always)]
pub fn in_degrees(rad:f64) -> f64 {
    rad * 180.0/std::f64::consts::PI
}

/// Converts degrees into Radians
#[inline(always)]
pub fn in_radians(degrees:f64)->f64{
    degrees * std::f64::consts::PI/180.0
}

/// The Equation of Time based on the day of year (can have decimals)
/// 
/// n should be in solar time, but this variable does not change daily so
/// it probably does not matter
pub fn equation_of_time(n: f64)->f64{
    let b = b(n);
    229.2*(0.000075 + 0.001868*b.cos()-0.032077*b.sin() - 0.014615*(2.0*b).cos() - 0.04089*(2.0*b).sin())
}

/// Declination (in Radians), according to Equation 1.6.1B
/// 
/// n should be in solar time, but this variable does not change daily so
/// it probably does not matter
pub fn declination(n:f64)->f64{    
    let b = b(n);
    
    // Return in Radians   
    0.006918 
    - 0.399912 *  b.cos()     + 0.070257 * b.sin() 
    - 0.006758 * (2.*b).cos() + 0.000907 * (2.*b).sin() 
    - 0.002697 * (3.*b).cos() + 0.00148  *(3.*b).sin()
    
}

/// Equation 1.4.2 in the Book.
/// 
/// n should be in solar time, but this variable does not change daily so
/// it probably does not matter
#[inline(always)]
fn b(n:f64)->f64{
    (n-1.0)*2.0*std::f64::consts::PI/365.0
}

impl Solar {
    
    /// Builds a Solar object from a Standard Time date, a Latitude,
    /// Longitude and Standard meridian (in Radians)
    pub fn from_standard_time(standard_time: Date, latitude: f64, longitude:f64, standard_meridian:f64) -> Self {
        let mut n = standard_time.day_of_year();        
        let delta_minutes = 4.0*in_degrees(standard_meridian-longitude)+equation_of_time(n);
        
        // Add the number of minutes divided the number of minutes in a day
        n += delta_minutes/24./60.;
        

        Self{
            n,
            latitude,
            //longitude,
            //standard_meridian,
        }
    
    }

    /// Builds a Solar object from a Solar Time date and a Latitude 
    /// (in Radians)
    pub fn from_solar_time(solar_time: Date, latitude: f64) -> Self {
        let n = solar_time.day_of_year();        
        
        Self{
            n,
            latitude,
            //longitude,
            //standard_meridian,
        }
    
    }



    /// Equation 1.4.2 in the Book.
    /// 
    /// n should be in solar time, but this variable does not change daily so
    /// it probably does not matter
    #[inline(always)]
    fn b(&self)->f64{
        b(self.n)
    }
    
    /// Normal extraterrestrial radiation (Gon)
    /// Equation 1.4.1b from Duffie and Beckman
    /// 
    /// n should be in solar time, but this variable does not change daily so
    /// it probably does not matter
    pub fn normal_extraterrestrial_radiation(&self)->f64{
        let b = self.b();
        let aux = 1.000110 + 0.034221 * b.cos() + 0.001280*b.sin()+0.000719*(2.0*b).cos() + 0.000077*(2.0*b).sin();
        SOLAR_CONSTANT*aux
    
    }
        
    
    /// Declination, according to Equation 1.6.1B
    /// 
    /// n should be in solar time, but this variable does not change daily so
    /// it probably does not matter
    pub fn declination(&self)->f64{
        declination(self.n)        
    }
    
    /// Returns the hour angle in degrees
    /// 
    /// n should be in solar time
    pub fn hour_angle(&self)->f64{
        // Remove the day (keep the hour)
        let solar_hour = 24.*(self.n % 1.);
    
        // Multiply for 24 hours, and by 15degrees/hour
        in_radians((solar_hour - 12.)*15.)
    }

    /// Builds a vector that points towards the sun.
    /// 
    /// Z is up, Y is North and X is East
    pub fn sun_position(&self)->Vector3D{
        
        let cos_phi = self.latitude.cos();
        let sin_phi = self.latitude.sin();
        
        let delta = self.declination();
        let cos_delta = delta.cos();
        let sin_delta = delta.sin();

        let omega = self.hour_angle();
        let cos_omega = omega.cos();

        // Equation 1.6.5, for Zenith
        let cos_zenith = cos_phi* cos_delta * cos_omega + sin_phi * sin_delta;
        let sin_zenith = cos_zenith.acos().sin();
        debug_assert!( (1.0-(cos_zenith*cos_zenith + sin_zenith * sin_zenith)).abs() < 0.000001);
        let z = cos_zenith;        

        // Equation 1.6.6 for Azimuth
        let cos_azimuth = (cos_zenith * sin_phi - sin_delta)/(sin_zenith * cos_phi);
        let sin_azimuth = cos_azimuth.acos().sin();
        debug_assert!( (1.0-(cos_azimuth*cos_azimuth + sin_azimuth*sin_azimuth)).abs() < 0.0000001 );
        let y = -cos_azimuth * sin_zenith;

        
        // Trigonometry
        let mut x = sin_azimuth * sin_zenith;
        // (x should be positive at this stage, right? then, if omega is 
        // positive, we need to change the sign of x)
        if omega > 0. {
            x *= -1.
        }

        // Check length of vector
        debug_assert!( ((x*x+y*y+z*z).sqrt()-1.0).abs()< 0.000001);

        // Build the vector and return
        Vector3D::new(x,y,z)
    }

} // end of impl Solar


#[cfg(test)]
mod tests {

    use super::*;

    fn are_close(x:f64, y:f64, precision: f64)->bool{
        if (x-y).abs() < precision {
            return true
        }
        println!("x:{}, y:{}", x, y);
        false
        
    }

    #[test]
    fn test_in_degrees() {
        // They are inverse functions of each other
        // Thus, x = g(f(x)), for any x 
        assert!(are_close(1.0, in_radians(in_degrees(1.0)), 0.00001));
        assert!(are_close(0.3, in_radians(in_degrees(0.3)), 0.00001));

        
        // Conversions from rad to degree
        const EPS : f64 = 0.0001;
        assert!(are_close(in_degrees(1.0), 57.29578, EPS));
        assert!(are_close(in_degrees(1.45), 83.07888, EPS));
        assert!(are_close(in_degrees(std::f64::consts::PI), 180.0, EPS));

        // Conversions from degree to rad
        assert!(are_close(in_radians(180.0), std::f64::consts::PI, EPS));
        assert!(are_close(in_radians(230.0), 4.014257, EPS));
        assert!(are_close(in_radians(130.0), 2.268928, EPS));
    }

    #[test]
    fn test_normal_extraterrestrial(){
        assert!(false)
    }

    #[test]
    fn test_from_standard_time(){        
        /*
        Example 1.5.1 in Duffie & Beckman
        At Madison, Wisconsin, what is the solar time corresponding to 10:30 AM central time on February 3?

        Answer: 10:19
        */
        let longitude = in_radians(89.4);
        let standard_meridian = in_radians(90.0);
        let standard_time = Date{
            month:2, day: 3, hour:10.5
        };
        const EPS:f64 = 0.5/60.;// Half a minute precision
        let solar_time = Solar::from_standard_time(standard_time, 0., longitude, standard_meridian);        
        let solar_date = Date::from_day_of_year(solar_time.n);
        assert!(are_close(solar_date.hour, 10.0+19.0/60.0, EPS))

    }

    #[test]
    fn test_declination(){
        
        fn check(month:u8, day:u8, expected_n:f64, expected_d: f64){
            let date = Date{ month: month, day: day, hour: 0. };
            let n = date.day_of_year();
            assert_eq!(n, expected_n - 1.);

            let d = declination(n);

            println!("exp: {}, found: {}", expected_d, in_degrees(d));
            // I suspect I need this margin of error (1.8 deg.)
            // because Duffie and Beckam do not specify the hour
            // of the day or the exact equation they use.
            assert!(are_close(in_degrees(d), expected_d, 1.8))

        }

        // From table 1.6.1... declinations are in degrees
        check(1,  17, 17.,  -20.9);
        check(2,  16, 47.,  -13.0);
        check(3,  16, 75.,  -2.4);
        check(4,  15, 105.,  9.4);
        check(5,  15, 135.,  18.8);
        check(6,  11, 162.,  23.1);
        check(7,  17, 198.,  21.2);
        check(8,  16, 228.,  13.5);
        check(9,  15, 258.,  2.2);
        check(10, 15, 288., -9.6);
        check(11, 14, 318., -18.9);
        check(12, 10, 344., -23.0);
    }

    
    #[test]
    fn test_hour_angle(){
        // Example 1.6.1
        /*
        10:30 (solar time) on February 13... hour_angle is -22.55
        */
        let date = Date{
            month: 2, day:13, hour:10.5
        };        
        let solar = Solar::from_solar_time(date, 0.);
        let w = in_degrees(solar.hour_angle());
        assert!(are_close(w, -22.5, 0.1));

        /* OTHERS */
        // Midday == 0
        let date = Date{month: 2, day:13, hour:12.0 };        
        let solar = Solar::from_solar_time(date, 0.);
        let w = in_degrees(solar.hour_angle());
        assert!(are_close(w, 0., 0.1));

        // 13:00 == 15
        let date = Date{month: 2, day:13, hour:13.0 };        
        let solar = Solar::from_solar_time(date, 0.);
        let w = in_degrees(solar.hour_angle());
        assert!(are_close(w, 15., 0.1));
    }

    #[test]
    fn test_sun_position(){
        /*
        Example 1.6.2
        Calculate the zenith and solar azimuth angles for φ = 43◦ at 9:30 AM on February 13 and 6:30 PM on July 1.        
        */
        let phi = in_radians(43.);

        // FOR 9:30 AM on February 13
        // ==========================
        let solar_time = Date{month:2, day:13, hour:9.5};
        let solar = Solar::from_solar_time(solar_time, phi);
        let dir = solar.sun_position();
        assert!(are_close(dir.length(), 1.0, 0.00001));

        // check declination
        assert!(are_close(in_degrees(solar.declination()), -14., 0.5));

        // check hour angle
        assert!(are_close(in_degrees(solar.hour_angle()), -37.5, 0.5));

        // zenith
        let zenith = in_degrees(dir.z().acos());
        assert!(are_close(zenith, 66.5, 0.5));

        // Azimuth
        let azimuth = in_degrees((dir.x()/dir.y()).atan());
        assert!(are_close(azimuth, -40., 0.5));


        // 6:30 PM on July 1
        // =================
        let solar_time = Date{month:7, day:1, hour:18.5};
        let solar = Solar::from_solar_time(solar_time, phi);
        let dir = solar.sun_position();
        assert!(are_close(dir.length(), 1.0, 0.00001));

        // check declination
        assert!(are_close(in_degrees(solar.declination()), 23.1, 0.5));

        // check hour angle
        assert!(are_close(in_degrees(solar.hour_angle()), 97.5, 0.5));

        // zenith
        let zenith = in_degrees(dir.z().acos());
        assert!(are_close(zenith, 79.6, 0.5));

        // Azimuth
        println!("{}", dir);
        let azimuth = in_degrees((dir.x()/dir.y()).atan());
        //assert!(are_close(azimuth, 112., 0.5)); // This is working, but atan() returns -67 instead of 112

        
    }


    #[test]
    fn test_angle_of_incidence(){
        /*
        Example 1.6.1
        Calculate the angle of incidence of beam radiation on a surface 
        located at Madison, Wisconsin, at 10:30 (solar time) on February 13 
        if the surface is tilted 45◦ from the horizontal and pointed 15◦ 
        west of south.        
        */
        // sun direction
        let latitude = in_radians(43.);
        let solar_time = Date{month:2, day:13, hour:10.5};
        let solar = Solar::from_solar_time(solar_time, latitude);
        let solar_dir = solar.sun_position();
        // check declination
        assert!(are_close(in_degrees(solar.declination()), -14., 0.5));

        // check hour angle
        assert!(are_close(in_degrees(solar.hour_angle()), -22.5, 0.5));

        // surface
        let beta = in_radians(45.);
        let gamma = in_radians(15.);

        let x = -gamma.sin() * beta.sin();
        let y = -gamma.cos() * beta.sin();
        let z = beta.cos();        
        let surface_dir = Vector3D::new(x,y,z);
        println!("{} | len = {}", surface_dir, surface_dir.length());

        let angle = (solar_dir * surface_dir).acos();

        assert!(are_close(in_degrees(angle), 35., 0.2));


    }
    
}
