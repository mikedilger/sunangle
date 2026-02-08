use astrotime::{DateTime, Duration, Epoch, Gregorian, Instant, Tt};

mod types;
pub use types::*;

mod vec;
pub use vec::*;

// NOTE:  f64 gives at least 15 digits of precision.
//        but many formulas don't offer that.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ObliquityFormula {
    // Used until 1983, Newcomb was used that
    // used B1900.0 and Tropical centuries. We
    // skip that one.
    /// Jet Propulsion Labs, 1984
    JPL1984,

    /// NOAA JPL1984 corrected
    /// From their spreadsheet
    NOAA,

    /// International Astronomical Union resolution in 2006
    /// to adopt PO3, high resolution over a few hundred years
    /// around the current time.
    #[default]
    IAU2006PO3,

    /// Laksar, intended to span longer time periods (over
    /// 1000 years).
    Laksar,
}

// We used data here:
// https://gml.noaa.gov/grad/solcalc/calcdetails.html
// which reference Astronomical Algorithms by Jean Meeus.
//
// Also, James Still's blog posts starting in Feb 2019
// reference the same source and explain things in much
// more detail:  https://squarewidget.com/julian-day/

/// Coordinates (when and where) for the calculations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coords {
    pub latitude: Degrees,
    pub longitude: Degrees,
    pub time: Instant,
    pub obliquity_formula: ObliquityFormula,
}

impl Coords {
    /// Julian day (days since the Julian epoch)
    #[inline]
    pub fn julian_day(&self) -> f64 {
        let duration: Duration = self.time - Epoch::JulianPeriod.as_instant();
        let seconds = duration.seconds_part() as f64
            + (duration.attos_part() / astrotime::ATTOS_PER_SEC_I64) as f64;
        seconds / 86400.0
    }

    /// Julian century, from the J2000.0 epoch
    #[inline]
    pub fn julian_century_j2000(&self) -> f64 {
        // This converts julian day to J2000.0,
        // then to julian centuries
        (self.julian_day() - 2451545.) / 36525.
    }

    /// This is the mean longitude of the sun, corrected for the abberation
    /// of light, in the ecliptic coordinate system.
    pub fn solar_geometric_mean_longitude(&self) -> Degrees {
        let jc = self.julian_century_j2000();
        Degrees((280.46646 + jc * (36000.76983 + jc * 0.0003032)).rem_euclid(360.))
    }

    /// This is the mean anomoly of the sun (pretending the sun orbits the
    /// earth) in the ecliptic coordinate system.
    /// The orbit is an ellipse, and the anomoly measures the angle swept out since
    /// the perihelion (the closest approach).
    pub fn solar_geometric_mean_anomoly(&self) -> Degrees {
        let jc = self.julian_century_j2000();
        Degrees(357.52911 + jc * (35999.05029 - 0.0001537 * jc))
    }

    /// The ecentricity of Earth's orbit about the sun (how flat the ellipse is
    /// compared to circular which would be 0.0)
    pub fn earth_orbit_eccentricity(&self) -> f64 {
        let jc = self.julian_century_j2000();
        0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)
    }

    pub fn sun_equation_of_center(&self) -> Degrees {
        let jc = self.julian_century_j2000();
        let ma = self.solar_geometric_mean_anomoly();

        Degrees(
            sin(ma) * (1.914602 - jc * (0.004817 + 0.000014 * jc))
                + sin(ma * 2.0) * (0.019993 - 0.000101 * jc)
                + sin(ma * 3.0) * 0.000289,
        )
    }

    /// This is the sun's true longitude in the ecliptic coordinate system.
    pub fn sun_true_longitude(&self) -> Degrees {
        self.solar_geometric_mean_longitude() + self.sun_equation_of_center()
    }

    /// This is the sun's true anomoly in the ecliptic coordinate system.
    /// The orbit is an ellipse, and the anomoly measures the angle swept out since
    /// the perihelion (the closest approach).
    pub fn sun_true_anomoly(&self) -> Degrees {
        self.solar_geometric_mean_anomoly() + self.sun_equation_of_center()
    }

    /// This is the radius vector magnitude, the distance from the center
    /// of the Sun to the center of the Earth, in AU.
    pub fn sun_rad_vector_au(&self) -> f64 {
        let ee = self.earth_orbit_eccentricity();
        let ta = self.sun_true_anomoly();
        (1.000001018 * (1.0 - ee * ee)) / (1.0 + ee * cos(ta))
    }

    /// This corrects the Sun's longitude for perturbations that affect
    /// it's apparent longitude
    pub fn sun_apparent_longitude(&self) -> Degrees {
        let jc = self.julian_century_j2000();
        let tl = self.sun_true_longitude();
        Degrees(tl.0 - 0.00569 - 0.00478 * sin(Degrees(125.04 - 1934.136 * jc)))
    }

    /// This is the axis tilt of the Earth with respect to the celestial
    /// equator
    pub fn obliquity(&self) -> Degrees {
        let jc = self.julian_century_j2000();

        match self.obliquity_formula {
            // Newcomb, pre 1984 (skipped)
            // ε = 23°27′8.26″ − 46.845″ T − 0.0059″ T2 + 0.00181″ T3
            ObliquityFormula::JPL1984 => {
                // ε = 23°26′21.448″ − 46.8150″ T − 0.00059″ T2 + 0.001813″ T3
                Degrees::from_dms(
                    23,
                    26,
                    21.448 - 46.8150 * jc - 0.00059 * jc.powi(2) + 0.001813 * jc.powi(3),
                )
            }
            ObliquityFormula::NOAA => {
                let jc = self.julian_century_j2000();
                Degrees::from_dms(
                    23,
                    26,
                    21.448 - 46.8150 * jc - 0.00059 * jc.powi(2) + 0.001813 * jc.powi(3),
                ) + Degrees(0.00256 * cos(Degrees(125.04 - 1934.136 * jc)))
            }
            ObliquityFormula::IAU2006PO3 => {
                // ε = 23°26′21.406″ − 46.836769″ T − 0.0001831″ T2 + 0.00200340″ T3 − 5.76″ × 10−7 T4 − 4.34″ × 10−8 T5
                Degrees::from_dms(
                    23,
                    26,
                    21.406 - 46.836769 * jc - 0.0001831 * jc.powi(2) + 0.00200340 * jc.powi(3)
                        - 0.000000576 * jc.powi(4)
                        - 0.0000000434 * jc.powi(5),
                )
            }
            ObliquityFormula::Laksar => {
                // ε = 23°26′21.448″ − 4680.93″ t − 1.55″ t2 + 1999.25″ t3 − 51.38″ t4 − 249.67″ t5 − 39.05″ t6 + 7.12″ t7 + 27.87″ t8 + 5.79″ t9 + 2.45″ t10
                let t = jc / 100.0;
                Degrees::from_dms(
                    23,
                    26,
                    21.448 - 4680.93 * t - 1.55 * t.powi(2) + 1999.25 * t.powi(3)
                        - 51.38 * t.powi(4)
                        - 249.67 * t.powi(5)
                        - 39.05 * t.powi(6)
                        + 7.12 * t.powi(7)
                        + 27.87 * t.powi(8)
                        + 5.79 * t.powi(9)
                        + 2.45 * t.powi(10),
                )
            }
        }
    }

    /// This is the right ascension of the Sun
    pub fn sun_right_ascension(&self) -> Degrees {
        let sal = self.sun_apparent_longitude();
        let oc = self.obliquity();
        atan2(cos(oc) * sin(sal), cos(sal))
    }

    /// This is the declination of the Sun
    pub fn sun_declination(&self) -> Degrees {
        let sal = self.sun_apparent_longitude();
        let oc = self.obliquity();
        asin(sin(oc) * sin(sal))
    }

    fn var_y(&self) -> f64 {
        let oc2 = self.obliquity() / 2.0;
        let toct = tan(oc2);
        toct * toct
    }

    pub fn equation_of_time(&self) -> f64 {
        let vary = self.var_y();
        let sgmlat = self.solar_geometric_mean_longitude();
        let sgmanom = self.solar_geometric_mean_anomoly();
        let orbitecc = self.earth_orbit_eccentricity();
        4.0 * (vary * (2.0 * sgmlat.0.to_radians()).sin() - 2.0 * orbitecc * sin(sgmanom)
            + 4.0 * orbitecc * vary * sin(sgmanom) * (2.0 * sgmlat.0.to_radians()).cos()
            - 0.5 * vary * vary * (4.0 * sgmlat.0.to_radians()).sin()
            - 1.25 * orbitecc * orbitecc * (2.0 * sgmanom.0.to_radians()).sin())
        .to_degrees()
    }

    pub fn ha_sunrise(&self) -> Degrees {
        let decl = self.sun_declination();
        acos(
            cos(Degrees(90.833)) / (cos(self.latitude) * cos(decl))
                - tan(self.latitude) * tan(decl),
        )
    }

    pub fn solar_noon_utc(&self) -> Hours {
        let eot = self.equation_of_time();
        let day_frac = (720.0 - 4.0 * self.longitude.0 - eot) / 1440.0;
        Hours(day_frac * 24.0)
    }

    pub fn sunrise_time_utc(&self) -> Hours {
        let has = self.ha_sunrise();
        let sn = self.solar_noon_utc();
        let mut h = Hours(sn.0 - 24.0 * (has.0 * 4.0) / 1440.0);
        h.normalize();
        h
    }

    pub fn sunset_time_utc(&self) -> Hours {
        let has = self.ha_sunrise();
        let sn = self.solar_noon_utc();
        let mut h = Hours(sn.0 + 24.0 * (has.0 * 4.0) / 1440.0);
        h.normalize();
        h
    }

    pub fn sunlight_duration_minutes(&self) -> f64 {
        let has = self.ha_sunrise();
        has.0 * 8.0
    }

    pub fn true_solar_time_minutes(&self) -> f64 {
        let eot = self.equation_of_time();
        let dt: DateTime<Gregorian, Tt> = self.time.into();
        let day_frac = dt.day_fraction();

        (day_frac * 1440.0 + eot + 4.0 * self.longitude.0).rem_euclid(1440.0)
    }

    pub fn hour_angle(&self) -> Degrees {
        let tst = self.true_solar_time_minutes();
        if tst / 4.0 < 0.0 {
            Degrees(tst / 4.0 + 180.0)
        } else {
            Degrees(tst / 4.0 - 180.0)
        }
    }

    pub fn solar_zenith_angle(&self) -> Degrees {
        let decl = self.sun_declination();
        let ha = self.hour_angle();
        acos(sin(self.latitude) * sin(decl) + cos(self.latitude) * cos(decl) * cos(ha))
    }

    pub fn solar_elevation_angle(&self) -> Degrees {
        Degrees(90.0) - self.solar_zenith_angle()
    }

    pub fn approx_atmospheric_refraction(&self) -> Degrees {
        let elev = self.solar_elevation_angle();
        let tanh = elev.0.to_radians().tan();

        Degrees(
            match elev.0 {
                x if x > 85.0 => 0.0,
                x if x > 5.0 => 58.1 / tanh - 0.07 / tanh.powi(3) + 0.000086 / tanh.powi(5),
                x if x > -0.575 => {
                    1735.0
                        + elev.0 * (-518.2 + elev.0 * (103.4 + elev.0 * (-12.79 + elev.0 * 0.711)))
                }
                _ => -20.772 / tanh,
            } / 3600.0,
        )
    }

    pub fn solar_elevation_corrected(&self) -> Degrees {
        self.solar_elevation_angle() + self.approx_atmospheric_refraction()
    }

    pub fn solar_zenith_corrected(&self) -> Degrees {
        self.solar_zenith_angle() - self.approx_atmospheric_refraction()
    }

    // Clockwise from N
    pub fn solar_azimuth_angle(&self) -> Degrees {
        let ha = self.hour_angle();
        let decl = self.sun_declination();
        let za = self.solar_zenith_angle();

        let d = if ha.0 > 0.0 {
            acos((sin(self.latitude) * cos(za) - sin(decl)) / (cos(self.latitude) * sin(za))).0
                + 180.0
        } else {
            540.0
                - acos((sin(self.latitude) * cos(za) - sin(decl)) / (cos(self.latitude) * sin(za)))
                    .0
        };

        Degrees(d.rem_euclid(360.0))
    }

    // A vector pointing towards where the solar irradiance appears to be
    // coming from.
    pub fn sun_vector(&self) -> Vec3 {
        let z = self.solar_zenith_corrected();
        let a = self.solar_azimuth_angle();

        Vec3::new(
            sin(z) * sin(a),
            sin(z) * cos(a),
            cos(z)
        )
    }
}

#[cfg(test)]
mod newtest {
    use crate::*;

    #[test]
    fn test_coords() {
        // This uses the first line data from the NOAA spreadsheet to verify
        // that our formulas generate very close to the same output (usually
        // 14 bits of precision), and the expected values were modified
        // to match the precise f64 data we are getting so we could use
        // equality testing between floats.

        let c = Coords {
            latitude: Degrees(40.0),
            longitude: Degrees(-105.0),
            time: DateTime::<Gregorian, Tt>::new(2010, 6, 21, 7, 6, 0, 0)
                .unwrap()
                .into(),
            obliquity_formula: ObliquityFormula::NOAA,
        };

        let julian_day = c.julian_day();
        println!("julian day {julian_day}");
        assert_eq!(julian_day, 2455368.7958333334);

        let julian_century = c.julian_century_j2000();
        println!("julian century (J2000.0) {julian_century}");
        assert_eq!(julian_century, 0.10468982432124285);

        let mean_longitude = c.solar_geometric_mean_longitude();
        println!("mean longitude {mean_longitude}");
        assert_eq!(mean_longitude.0, 89.38073225525932);

        let mean_anomoly = c.solar_geometric_mean_anomoly();
        println!("mean anomoly {mean_anomoly}");
        assert_eq!(mean_anomoly.0, 4126.263358907141);

        let earth_orbit_eccentricity = c.earth_orbit_eccentricity();
        println!("earth orbit eccentricity {earth_orbit_eccentricity}");
        assert_eq!(earth_orbit_eccentricity, 0.016704231765228162);

        let sun_eq_of_ctr = c.sun_equation_of_center();
        println!("sun equation of center {sun_eq_of_ctr}");
        assert_eq!(sun_eq_of_ctr.0, 0.44549228519368134);

        let sun_true_long = c.sun_true_longitude();
        println!("sun true longitude {sun_true_long}");
        assert_eq!(sun_true_long.0, 89.826224540453);

        let sun_true_anom = c.sun_true_anomoly();
        println!("sun true anomoly {sun_true_anom}");
        assert_eq!(sun_true_anom.0, 4126.708851192335);

        let sun_rad_vec = c.sun_rad_vector_au();
        println!("sun rad vec (AU) {sun_rad_vec}");
        assert_eq!(sun_rad_vec, 1.0162428418623726);

        let sun_app_long = c.sun_apparent_longitude();
        println!("sun apparent longitude {sun_app_long}");
        assert_eq!(sun_app_long.0, 89.82520022844825);

        let obliquity = c.obliquity();
        println!("obliquity {obliquity}");
        assert_eq!(obliquity.0, 23.438486218294376);

        let sun_right_ascension = c.sun_right_ascension();
        println!("sun right ascension {sun_right_ascension}");
        assert_eq!(sun_right_ascension.0, 89.80948008441563);

        let sun_declination = c.sun_declination();
        println!("sun declination {sun_declination}");
        assert_eq!(sun_declination.0, 23.438370619286854);

        let var_y = c.var_y();
        println!("var_y {var_y}");
        assert_eq!(var_y, 0.043031489687856965);

        let eot = c.equation_of_time();
        println!("equation of time {eot}");
        assert_eq!(eot, -1.7153675784683604);

        let ha_sunrise = c.ha_sunrise();
        println!("HA sunrise {ha_sunrise})");
        assert_eq!(ha_sunrise.0, 112.61041006988431);

        let solar_noon_utc = c.solar_noon_utc();
        println!("solar noon UTC {solar_noon_utc}");
        assert_eq!(solar_noon_utc.0, 19.02858945964114);

        let sunrise = c.sunrise_time_utc();
        println!("sunrise UTC {sunrise}");
        assert_eq!(sunrise.0, 11.521228788315518); // 11h 31m 16.42363  (4:31:16 at -7)

        let sunset = c.sunset_time_utc();
        println!("sunset UTC {sunset}");
        assert_eq!(sunset.0, 2.5359501309667607); // 2h 32m 9.42047  (19:32:09 at -7)

        let sunlight_minutes = c.sunlight_duration_minutes();
        println!("sunlight (min) {sunlight_minutes}");
        assert_eq!(sunlight_minutes, 900.8832805590745);

        let true_solar_time_min = c.true_solar_time_minutes();
        println!("True solar time (TT, minutes) {true_solar_time_min}");
        assert_eq!(true_solar_time_min, 4.284632421531626);

        let hour_angle = c.hour_angle();
        println!("Hour angle {hour_angle}");
        assert_eq!(hour_angle.0, -178.9288418946171);

        let solar_zenith_angle = c.solar_zenith_angle();
        println!("Solar zenith angle {solar_zenith_angle}");
        assert_eq!(solar_zenith_angle.0, 116.55376211918129);

        let solar_elevation_angle = c.solar_elevation_angle();
        println!("Solar elevation angle {solar_elevation_angle}");
        assert_eq!(solar_elevation_angle.0, -26.553762119181286);

        let approx_refr = c.approx_atmospheric_refraction();
        println!("Approx atmos refraction {approx_refr}");
        assert_eq!(approx_refr.0, 0.01154568659192415);

        let solar_elevation_corrected = c.solar_elevation_corrected();
        println!("Solar elevation corrected {solar_elevation_corrected}");
        assert_eq!(solar_elevation_corrected.0, -26.542216432589363);

        let solar_azimuth_angle = c.solar_azimuth_angle();
        println!("Solar azimuth angle {solar_azimuth_angle}");
        assert_eq!(solar_azimuth_angle.0, 1.0986711838851306);
    }
}

// Astronomers tend to use degrees, but computer trig functions use radians
// so we wrap those functions
fn sin(d: Degrees) -> f64 {
    d.0.to_radians().sin()
}
fn cos(d: Degrees) -> f64 {
    d.0.to_radians().cos()
}
fn tan(d: Degrees) -> f64 {
    d.0.to_radians().tan()
}
fn asin(x: f64) -> Degrees {
    Degrees(x.asin().to_degrees())
}
fn acos(x: f64) -> Degrees {
    Degrees(x.acos().to_degrees())
}
fn atan2(x: f64, y: f64) -> Degrees {
    Degrees(x.atan2(y).to_degrees())
}

/* -----------------------

BELOW IS OLDER CODE, but might be newer formulas.
It needs to be reviewed.


fn julian_century_j2000(instant: Instant) -> f64 {
    // Without astrotime we could compute this as follows:
    // spreadsheet_date_epoch = 2415018.5;
    // julian_day = spreadsheet_date + spreadsheet_date_epoch
    //              + hours/24 - timezoneoffsethours/24

    let duration: Duration = instant - Epoch::J2000_0.as_instant();
    let seconds = duration.seconds_part() as f64;
    seconds / 3155760000.0  // (86400.0 * 36525.0)
}


/// The solar altitude and azimuth angles at the given time and
/// location.
pub fn sun_altitude_azimuth(
    instant: Instant,
    latitude: Degrees,
    longitude: Degrees
) -> (Degrees, Degrees) {
    let ephemeris = sun_ephemeris(instant);
    let gst = greenwich_sidereal_time(instant);

    let ra = ephemeris.right_ascension.0;
    let d = ephemeris.declination;

    let lha = Degrees(gst + longitude.0 - ra);

    let alt = Degrees(sin(latitude) * sin(d)
        + cos(latitude) * cos(d) * cos(lha));

    let azimuth = atan2(
        -sin(lha),
        tan(d)*cos(latitude) - sin(latitude)*cos(lha)
    );

    (alt, azimuth)
}


// https://lweb.cfa.harvard.edu/~jzhao/times.html
fn greenwich_sidereal_time(instant: Instant) -> f64 {
    let t = julian_century_j2000(instant);

    24110.54841
        + 8640184.812866 * t
        + 0.093104 * t*t
        - 0.0000062 * t*t*t
}

/// This specifies the angle of an astronomical object
/// relative to the Earth, and is not specific to any
/// location on Earth.
pub struct Ephemeris {
    // Degrees north or south of the celestial equator
    pub right_ascension: Hours,
    pub declination: Degrees,
}

impl fmt::Display for Ephemeris {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RA: {}, Decl: {}",
            self.right_ascension, self.declination
        )
    }
}

/// Get the apparent position of the Sun relative to Earth
/// at some instant in time.
///
/// See https://aa.usno.navy.mil/faq/sun_approx
pub fn sun_ephemeris(instant: Instant) -> Ephemeris {
    let jcentury = julian_century_j2000(instant);
    println!("julian century = {jcentury}");

    let mean_longitude = (280.46646 + jcentury * (36000.76983 + jcentury * 0.0003032)).rem_euclid(360.0);

    println!("Mean lon = {mean_longitude}");

    let mean_anomoly = 357.52911 + jcentury * (35999.05029 - 0.0001537 * jcentury);

    println!("Mean anom = {mean_anomoly}");

    // This is adjusted for abberation
    let ecliptic_longitude = {
        let l = mean_longitude;
        let g = Degrees(mean_anomoly);
        Degrees(l + 1.915 * sin(g) + 0.020 * sin(g * 2.0))
    };

    // let beta = 0;
    // distance in AU
    // let r = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2.0*g);

    let obliquity = obliquity_hq(instant);

    let right_ascension = {
        let mut a = atan2(
            cos(obliquity) * sin(ecliptic_longitude),
            cos(ecliptic_longitude),
        );
        a.normalize();
        a
    };

    let declination = asin(sin(obliquity) * sin(ecliptic_longitude));

    Ephemeris {
        right_ascension: right_ascension.into(),
        declination: declination,
    }
}

// Axial tilt of the earth
fn obliquity_hq(instant: Instant) -> Degrees {
    // Simple formula
    // (23.439 - 0.0000004 * n).to_radians()

    let days_since_j2000 = {
        let duration: Duration = instant - Epoch::J2000_0.as_instant();
        duration.seconds_part() as f64 / 86400.0
    };

    // Laksar computation
    // julian 10 millenia
    let t = days_since_j2000 / 3652500.;

    // ε = 23° 26′ 21.448″ − 4680.93″ t − 1.55″ t2 + 1999.25″ t3 − 51.38″ t4 − 249.67″ t5 − 39.05″ t6 + 7.12″ t7 + 27.87″ t8 + 5.79″ t9 + 2.45″ t10    (23 + 26./60. + 21.448/3600.)
    Degrees(
        (23. + 26. / 60. + 21.448 / 3600.) - (4680.93 / 3600.) * t - (1.55 / 3600.) * t.powi(2)
            + (1999.25 / 3600.) * t.powi(3)
            - (51.38 / 3600.) * t.powi(4)
            - (249.67 / 3600.) * t.powi(5)
            - (39.05 / 3600.) * t.powi(6)
            + (7.12 / 3600.) * t.powi(7)
            + (27.87 / 3600.) * t.powi(8)
            + (5.79 / 3600.) * t.powi(9)
            + (2.45 / 3600.) * t.powi(10)
    )
}

// Axial tilt of the earth
// Used by most sun position calculators, but is less accurate than the above.
#[allow(dead_code)]
fn obliquity_approx(instant: Instant) -> Degrees {
    let days_since_j2000 = {
        let duration: Duration = instant - Epoch::J2000_0.as_instant();
        duration.seconds_part() as f64 / 86400.0
    };

    Degrees(23.439 - 0.00000036 * days_since_j2000)
}

#[cfg(test)]
mod test {
    use crate::{Degrees, Hours};
    use astrotime::{DateTime, Gregorian, Instant, Utc};

    #[test]
    fn test_sun_ephemeris() {
        struct Test {
            utc: Instant,
            right_ascension: Hours,
            declination: Degrees,
        }

        // http://www.ephemeris.com/ephemeris.php
        // Date/Time: 2026.01.28 06:33:22 UTC (GMT - Delta T), JD = 2461068.773171
        // Planet           Longitude       Latitude  Right Asc.  Declination
        // Sun            08 Aqr 20'28"     0°00'01"   20:43:03   -18°10'41"

        let test1 = Test {
            utc: DateTime::<Gregorian, Utc>::new(2026, 1, 28, 6, 33, 22, 0)
                .unwrap()
                .into(),
            // JD = 2461068.773171
            right_ascension: Hours::from_hms(20, 43, 03.),
            declination: Degrees::from_dms(-18, 10, 41.),
        };

        let sun_ephemeris = crate::sun_ephemeris(test1.utc);
        println!("{}", sun_ephemeris);
        // I'm getting accurate RA:
        // RA: 310.7631019629545° = 20h 43m 3.1444711090895794s,
        // But the decl is off by 6 seconds:
        // Decl: -18.17638664457702° = -18° -10’ -34.99192047726609”
    }

    #[test]
    fn test_sun_altitude_azimuth() {

        let latitude = Degrees(40);
        let longitude = Degrees(-105);
        let time: Instant =
            DateTime::<Gregorian, Utc>::new(2010, 6, 21, 0, 0, 0, 0)
            .unwrap()
            .into();


        // julian century = 0.26074663473
        //    I get         0.26074665659
        //    which is 69 seconds different.

        // mean long = 307.546061
        // mean anom = 9744.160316

        // azimuth 255.76546
        // altitude 11.50910

        let (altitude, azimuth) = crate::sun_altitude_azimuth(
            time,
            latitude,
            longitude
        );

        println!("Altitude: {altitude}");
        println!("Azimuth: {azimuth}");
}
}
*/
