/*
MIT License
Copyright (c) 2021 Germán Molina
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

const WHTEFFICACY: Float = 179.;
use crate::*;
use crate::{Float, PI};
use calendar::Date;
use geometry3d::Vector3D;
use matrix::Matrix;
use weather::CurrentWeather;

/// Specifies which units do we want returned
/// from the sky model
pub enum SkyUnits {
    /// Solar Radiance (i.e., W/m2.sr)
    Solar,

    /// Visible Radiance (i.e., W/m2.sr, but only visible)
    Visible,
}

/// This struct packs functions from two sources:
/// 1. Perez, R., Ineichen, P., Seals, R., Michalsky, J. and Stewart, R. (1990), "Modeling daylight availability and irradiance components from direct and global irradiance"
/// 2. Perez, R., R. Seals, and J. Michalsky (1993) All-Weather Model for Sky Luminance Distribution - Preliminary Configuration and Validation
/// I do not have access to the latter, so I borrowed code from Radiance's [`gendaymtx.c`](https://github.com/NREL/Radiance/blob/master/src/gen/gendaymtx.c) file
pub struct PerezSky {}

impl PerezSky {
    /// Equation 1
    fn sky_clearness(
        diffuse_horizontal_irrad: Float,
        direct_normal_irrad: Float,
        solar_zenith: Float,
    ) -> Float {
        const K: Float = 1.041; // Next to equation 1 of the paper
        let z_cubed = solar_zenith * solar_zenith * solar_zenith;
        ((diffuse_horizontal_irrad + direct_normal_irrad) / diffuse_horizontal_irrad + K * z_cubed)
            / (1. + K * z_cubed)
    }

    /// Equation 2
    fn sky_brightness(
        diffuse_horizontal_irrad: Float,
        air_mass: Float,
        extraterrestrial_irradiance: Float,
    ) -> Float {
        diffuse_horizontal_irrad * air_mass / extraterrestrial_irradiance
    }

    /// Equation 3... dew_point_temp in C
    fn precipitable_water_content(dew_point_temp: Float) -> Float {
        (0.07 * dew_point_temp - 0.075).exp()
    }

    /// Get the sky clearness category (Table 1 in the 1990 paper)
    fn clearness_category(clearness_index: Float) -> usize {
        if clearness_index < 1. {
            panic!(
                "A clearness index of {} is too low for Perez's Sky",
                clearness_index
            );
        } else if clearness_index < 1.065 {
            0
        } else if clearness_index < 1.230 {
            1
        } else if clearness_index < 1.5 {
            2
        } else if clearness_index < 1.95 {
            3
        } else if clearness_index < 2.8 {
            4
        } else if clearness_index < 4.5 {
            5
        } else if clearness_index < 6.2 {
            6
        } else if clearness_index < 12.01 {
            7
        } else {
            panic!(
                "A clearness index of {} is too high for Perez's Sky",
                clearness_index
            );
        }
    }

    // /// Calculates the direct illuminance/diffuse radiance ratio
    // /// according to Equation 8 and Table 4 of Perez et al. 1990
    fn direct_illuminance_ratio(p_water_content: Float, zenit: Float, sky_brightness: Float, index: usize)->Float{
        if index > 7 {
            panic!("Table 4 of Perez's paper has only 8 sky clearness categories (starting from 0)... received {}", index)
        }
        const TABLE : [(Float, Float, Float, Float); 8] =
        [
            (  57.20, -4.55, -2.98, 117.12 ),
            (  98.99, -3.46, -1.21,  12.38 ),
            ( 109.83, -4.90, -1.71,  -8.81 ),
            ( 110.34, -5.84, -1.99,  -4.56 ),
            ( 106.36, -3.97, -1.75,  -6.16 ),
            ( 107.19, -1.25, -1.51, -26.73 ),
            ( 105.75,  0.77, -1.26, -34.44 ),
            ( 101.18,  1.58, -1.10,  -8.29 )
        ];

        let v = TABLE[index].0 +
        TABLE[index].1*p_water_content +
        TABLE[index].2*(5.73*zenit - 5.).exp() +
        TABLE[index].3* sky_brightness;

        // return
        v.clamp(0., 9e9)
    }

    /// Calculates the diffuse illuminance/diffuse radiance ratio
    /// according to Equation 7 and Table 4 of Perez et al. 1990
    fn diffuse_illuminance_ratio(
        p_water_content: Float,
        cos_zenit: Float,
        sky_brightness: Float,
        index: usize,
    ) -> Float {
        if index > 7 {
            panic!("Table 4 of Perez's paper has only 8 sky clearness categories (starting from 0)... received {}", index)
        }

        const TABLE: [(Float, Float, Float, Float); 8] = [
            (97.24, -0.46, 12.00, -8.91),
            (107.22, 1.15, 0.59, -3.95),
            (104.97, 2.96, -5.53, -8.77),
            (102.39, 5.59, -13.95, -13.90),
            (100.71, 5.94, -22.75, -23.74),
            (106.42, 3.83, -36.15, -28.83),
            (141.88, 1.90, -53.24, -14.03),
            (152.23, 0.35, -45.27, -7.98),
        ];

        // Return
        TABLE[index].0
            + TABLE[index].1 * p_water_content
            + TABLE[index].2 * cos_zenit
            + TABLE[index].3 * sky_brightness.ln()
    }

    /// Calculates the Perez sky params
    fn calc_params(zenith: Float, epsilon: Float, mut delta: Float) -> [Float; 5] {
        const TABLE: [[Float; 20]; 8] = [
            /* Sky clearness (epsilon): 1.000 to 1.065 */
            [
                1.3525, -0.2576, -0.2690, -1.4366, -0.7670, 0.0007, 1.2734, -0.1233, 2.8000,
                0.6004, 1.2375, 1.0000, 1.8734, 0.6297, 0.9738, 0.2809, 0.0356, -0.1246, -0.5718,
                0.9938,
            ],
            /* Sky clearness (epsilon): 1.065 to 1.230 */
            [
                -1.2219, -0.7730, 1.4148, 1.1016, -0.2054, 0.0367, -3.9128, 0.9156, 6.9750, 0.1774,
                6.4477, -0.1239, -1.5798, -0.5081, -1.7812, 0.1080, 0.2624, 0.0672, -0.2190,
                -0.4285,
            ],
            /* Sky clearness (epsilon): 1.230 to 1.500 */
            [
                -1.1000, -0.2515, 0.8952, 0.0156, 0.2782, -0.1812, -4.5000, 1.1766, 24.7219,
                -13.0812, -37.7000, 34.8438, -5.0000, 1.5218, 3.9229, -2.6204, -0.0156, 0.1597,
                0.4199, -0.5562,
            ],
            /* Sky clearness (epsilon): 1.500 to 1.950 */
            [
                -0.5484, -0.6654, -0.2672, 0.7117, 0.7234, -0.6219, -5.6812, 2.6297, 33.3389,
                -18.3000, -62.2500, 52.0781, -3.5000, 0.0016, 1.1477, 0.1062, 0.4659, -0.3296,
                -0.0876, -0.0329,
            ],
            /* Sky clearness (epsilon): 1.950 to 2.800 */
            [
                -0.6000, -0.3566, -2.5000, 2.3250, 0.2937, 0.0496, -5.6812, 1.8415, 21.0000,
                -4.7656, -21.5906, 7.2492, -3.5000, -0.1554, 1.4062, 0.3988, 0.0032, 0.0766,
                -0.0656, -0.1294,
            ],
            /* Sky clearness (epsilon): 2.800 to 4.500 */
            [
                -1.0156, -0.3670, 1.0078, 1.4051, 0.2875, -0.5328, -3.8500, 3.3750, 14.0000,
                -0.9999, -7.1406, 7.5469, -3.4000, -0.1078, -1.0750, 1.5702, -0.0672, 0.4016,
                0.3017, -0.4844,
            ],
            /* Sky clearness (epsilon): 4.500 to 6.200 */
            [
                -1.0000, 0.0211, 0.5025, -0.5119, -0.3000, 0.1922, 0.7023, -1.6317, 19.0000,
                -5.0000, 1.2438, -1.9094, -4.0000, 0.0250, 0.3844, 0.2656, 1.0468, -0.3788,
                -2.4517, 1.4656,
            ],
            /* Sky clearness (epsilon): 6.200 to ... */
            [
                -1.0500, 0.0289, 0.4260, 0.3590, -0.3250, 0.1156, 0.7781, 0.0025, 31.0625,
                -14.5000, -46.1148, 55.3750, -7.2312, 0.4050, 13.3500, 0.6234, 1.5000, -0.6426,
                1.8564, 0.5636,
            ],
        ];

        let index = Self::clearness_category(epsilon);

        if epsilon > 1.065 && epsilon < 2.8 {
            if delta < 0.2 {
                delta = 0.2;
            }
        }

        /* Get Perez coefficients */
        let mut x = [[0.0; 4]; 5];
        for i in 0..5 {
            for j in 0..4 {
                x[i][j] = TABLE[index][4 * i + j];
            }
        }

        let mut perez_param = [0.0; 5];

        if index != 0 {
            /* Calculate parameter a, b, c, d and e (Eqn. 6) */
            for i in 0..5 {
                perez_param[i] = x[i][0] + x[i][1] * zenith + delta * (x[i][2] + x[i][3] * zenith);
            }
        } else {
            /* Parameters a, b and e (Eqn. 6) */
            perez_param[0] = x[0][0] + x[0][1] * zenith + delta * (x[0][2] + x[0][3] * zenith);
            perez_param[1] = x[1][0] + x[1][1] * zenith + delta * (x[1][2] + x[1][3] * zenith);
            perez_param[4] = x[4][0] + x[4][1] * zenith + delta * (x[4][2] + x[4][3] * zenith);

            /* Parameter c (Eqn. 7) */
            perez_param[2] = ((delta * (x[2][0] + x[2][1] * zenith)).powf(x[2][2])).exp() - x[2][3];

            /* Parameter d (Eqn. 8) */
            perez_param[3] =
                -(delta * (x[3][0] + x[3][1] * zenith)).exp() + x[3][2] + delta * x[3][3];
        }

        perez_param
    }

    /// Returns a function that can be used to calculate the Perez's sky
    /// values. The input to the returned function is a `Vector3D` pointing towards the
    /// sky (IT MUST BE NORMALIZED); the returned value is specified by the units parameter.    
    /// The Date must be in local time
    pub fn get_sky_func_standard_time(
        units: SkyUnits,
        solar: &Solar,
        date: Date,
        dew_point: Float,
        diffuse_horizontal_irrad: Float,
        direct_normal_irrad: Float,
        albedo: Float
    ) -> Box<dyn Fn(Vector3D) -> Float> {        

        // Convert local into solar time
        let day = Time::Standard(date.day_of_year());
        let sun_position = solar.sun_position(day).unwrap();

        // Convert day into solar time
        let day = solar.unwrap_solar_time(day);

        debug_assert!((1. - sun_position.length()) < 1e-5);

        let cos_zenit = sun_position.z;
        // If Dir is normalized, then its Z = cos(Zenith)
        let zenith = if cos_zenit <= 0. {
            // Limit zenith to 90 degrees
            PI / 2.
        } else if cos_zenit >= 0.9986295347545738 {
            // Limit Zenith to 3 degrees minimum
            /*
                The threshold above is equal to (3.*PI/180.).cos()
                would that have been optimized by the compiler?? I guess, but
                it did not allow me to create a constant of that value... so I
                did this just in case
            */
            (3. * PI / 180.).acos()
        } else {
            cos_zenit.acos()
        };

        let apwc = Self::precipitable_water_content(dew_point);
        let air_mass = air_mass(zenith);
        let extraterrestrial_irradiance = solar.normal_extraterrestrial_radiation(day);
        // RADIANCE SETS THESE LIMITS... I don't know if they are in the paper
        let sky_brightness = Self::sky_brightness(
            diffuse_horizontal_irrad,
            air_mass,
            extraterrestrial_irradiance,
        )
        .clamp(0.01, 9e9);
        let sky_clearness =
            Self::sky_clearness(diffuse_horizontal_irrad, direct_normal_irrad, zenith)
                .clamp(-9e9, 11.9);

        let index = Self::clearness_category(sky_clearness);

        let mut diff_illum = diffuse_horizontal_irrad
            * Self::diffuse_illuminance_ratio(apwc, cos_zenit, sky_brightness, index);
        let mut dir_illum = direct_normal_irrad * Self::direct_illuminance_ratio(apwc, zenith, sky_brightness, index);

        if let SkyUnits::Solar = units {
            diff_illum = diffuse_horizontal_irrad * WHTEFFICACY;
            dir_illum = direct_normal_irrad * WHTEFFICACY;
        }
        
        

        // Calculate Perez params
        let params = Self::calc_params(zenith, sky_clearness, sky_brightness);

        // Build closure
        let ret = move |dir: Vector3D| -> Float {
            debug_assert!((1. - dir.length()) < 1e-5);            

            let cosgamma = sun_position * dir;
            let gamma = cosgamma.acos();

            const MIN_DZ: Float = 0.01;

            let cos_zeta = if dir.z < MIN_DZ { MIN_DZ } else { dir.z };

            // return without the norm_diff_illum, because we need this closure
            // to calculate it.
            (1. + params[0] * (params[1] / cos_zeta).exp())
                * (1. + params[2] * (params[3] * gamma).exp() + params[4] * cosgamma * cosgamma)
        };

        // Calculate normalization coefficient
        let r = ReinhartSky::new(1);
        let mut norm_diff_illum = 0.;
        for i in 1..r.n_bins {
            let dir = r.bin_dir(i);

            let bin_dir = r.bin_dir(i);
            debug_assert!((1. - bin_dir.length()).abs() < 1e-5);
            norm_diff_illum += ret(dir) * r.bin_solid_angle(i) * bin_dir.z;
        }

        let norm_diff_illum = diff_illum / (norm_diff_illum * WHTEFFICACY);

        // Return
        Box::new(move |dir: Vector3D| -> Float { 
            // If ground
            if dir.z <= 0.0 {
                let mut global_horizontal = diff_illum;
                if cos_zenit > 0. {                    
                    global_horizontal += dir_illum * cos_zenit;
                }
                let ground = albedo * global_horizontal / PI/WHTEFFICACY;                
                return ground;
            }
            // If sky
            ret(dir) * norm_diff_illum 
        })
    }

    pub fn gen_sky_vec(
        mf: usize,
        solar: &Solar,
        date: Date,
        weather_data: CurrentWeather,
        units: SkyUnits,
        albedo: Float,
        add_sky: bool,
        add_sun: bool,
    ) -> Result<Matrix, String> {
        let r = ReinhartSky::new(mf);
        let mut vec = Matrix::new(0.0, r.n_bins, 1);
        Self::update_sky_vec(&mut vec, mf, solar, date, weather_data, units, albedo, add_sky, add_sun)?;
        Ok(vec)
    }

    pub fn update_sky_vec(
        vec: &mut Matrix,
        mf: usize,
        solar: &Solar,
        date: Date,
        weather_data: CurrentWeather,
        units: SkyUnits,
        albedo: Float,
        add_sky: bool,
        add_sun: bool,
    ) -> Result<(), String> {        

        let r = ReinhartSky::new(mf);
        let (rows, cols) = vec.size();
        debug_assert_eq!(cols, 1);
        if rows != r.n_bins {
            return Err(format!("when update_sky_vec() : number of elements of input vector ({}) does not match the number of bins for the Reinhart subdivition (MF {} require {} bins)", rows, mf, r.n_bins));
        }

        let dew_point = match weather_data.dew_point_temperature {
            None => {
                return Err(format!(
                    "Weather data needs dew point temperature for calculating sky-vec"
                ))
            }
            Some(v) => v,
        };
        let diffuse_horizontal_irrad = match weather_data.diffuse_horizontal_radiation {
            None => {
                return Err(format!(
                    "Weather data needs diffuse_horizontal_radiation for calculating sky-vec"
                ))
            }
            Some(v) => v,
        };
        let direct_normal_irrad = match weather_data.direct_normal_radiation {
            None => {
                return Err(format!(
                    "Weather data needs direct_normal_radiation for calculating sky-vec"
                ))
            }
            Some(v) => v,
        };

        // If it is nighttime, just fill with zeroes
        if direct_normal_irrad + diffuse_horizontal_irrad < 1e-4 {            
            for bin in 0..r.n_bins {                                
                vec.set(bin, 0, 0.0)?;
            }
            return Ok(())
        }

        if add_sky {
            let sky_func = Self::get_sky_func_standard_time(
                units,
                &solar,
                date,
                dew_point,
                diffuse_horizontal_irrad,
                direct_normal_irrad,
                albedo,
            );

            for bin in 0..r.n_bins {
                
                let dir = r.bin_dir(bin);
                let v = sky_func(dir);                
                vec.set(bin, 0, v)?;
            }
        }
        // Return
        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_gen_sky_vec() {
        let mf = 1;
        let lat = -41.41 * PI / 180.;        
        let lon = -174.87 * PI / 180.;
        let std_mer = -180. * PI / 180.;
        let month = 1;
        let day = 1;
        let hour = 5.5;
        let date = Date { month, day, hour };
        let solar = Solar::new(lat, lon, std_mer);
        let albedo = 0.2;
        let add_sky = true;
        let add_sun = false;

        let weather_data = CurrentWeather {
            dew_point_temperature: Some(11.),
            direct_normal_radiation: Some(538.),
            diffuse_horizontal_radiation: Some(25.),

            ..CurrentWeather::default()
        };
        let units = SkyUnits::Visible;

        let vec = PerezSky::gen_sky_vec(mf, &solar, date, weather_data, units, albedo, add_sky, add_sun).unwrap();
        println!("{}", vec);
    }

    #[test]
    fn test_visible_sky_func() {
        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 62.92923;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 55.68588;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 44.89257;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 34.8878;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 41.55878;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 38.9621;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 11.07143;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 37.18113;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 132.4133;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 37.42991;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 180.1197;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Visible,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 0 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 44.79652;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);
    }

    #[test]
    fn test_solar_sky_func() {
        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 83.71208;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 74.07656;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            0.017453292519943295,
            0.017453292519943295,
            0.017453292519943295,
        );
        let date = Date {
            month: 2,
            day: 2,
            hour: 12.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 500.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 2 2 12.0 -W 500.0 200.0 -a 1 -o 1 -m 1 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 59.71868;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 46.70042;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 55.63012;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(
            -0.5759586531581288,
            -0.6981317007977318,
            -0.6981317007977318,
        );
        let date = Date {
            month: 4,
            day: 5,
            hour: 10.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 600.0;
        let diffuse_horizontal_irrad = 200.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 4 5 10.0 -W 600.0 200.0 -a -33 -o -40 -m -40 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 52.15423;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 13.5277;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 45.43003;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(-0.8203047484373349, 0.8203047484373349, 0.8203047484373349);
        let date = Date {
            month: 11,
            day: 7,
            hour: 16.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 900.0;
        let diffuse_horizontal_irrad = 100.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 11 7 16.0 -W 900.0 100.0 -a -47 -o 47 -m 47 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 161.7902;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(
            0.04327423224079154,
            0.08654846448158308,
            0.9953073415382055,
        ));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.04327423224079154;  Dy=0.08654846448158308; Dz=0.9953073415382055' -f ./perezlum.cal -o '${intersky}'
        let expected = 58.99378;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(0.0, -0.9987523388778446, 0.04993761694389223));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=0.0;  Dy=-0.9987523388778446; Dz=0.04993761694389223' -f ./perezlum.cal -o '${intersky}'
        let expected = 283.889;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);

        let solar = Solar::new(0.8203047484373349, 0.20943951023931953, 0.20943951023931953);
        let date = Date {
            month: 1,
            day: 1,
            hour: 13.0,
        };
        let dew_point = 11.0;
        let direct_normal_irrad = 300.0;
        let diffuse_horizontal_irrad = 300.0;
        let fun = PerezSky::get_sky_func_standard_time(
            SkyUnits::Solar,
            &solar,
            date,
            dew_point,
            diffuse_horizontal_irrad,
            direct_normal_irrad,
            0.2,
        );
        let found = fun(Vector3D::new(-0.9912279006826347, 0.0, 0.13216372009101796));
        // Automatically generated using command: gendaylit 1 1 13.0 -W 300.0 300.0 -a 47 -o 12 -m 12 -O 1 | tail -n 1 | rcalc  -e 'A1=$2; A2=$3; A3=$4; A4=$5; A5=$6; A6=$7; A7=$8; A8=$9; A9=$10; A10=$11; Dx=-0.9912279006826347;  Dy=0.0; Dz=0.13216372009101796' -f ./perezlum.cal -o '${intersky}'
        let expected = 70.60438;
        println!(
            "{}, {}, {} ({}%)",
            expected,
            found,
            (expected - found).abs(),
            100. * (expected - found).abs() / expected
        );
        assert!(100. * (expected - found).abs() / expected < 3.);
    }
}
