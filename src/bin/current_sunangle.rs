use astrotime::*;
use sunangle::*;

fn main() {
    println!("Edit src/bin/sunangle.rs to change input parameters");

    let c = Coords {
        latitude: Degrees(0.0),
        longitude: Degrees(0.0),
        time: {
            let now = std::time::SystemTime::now();
            let now: Instant = TryFrom::try_from(now).unwrap();
            now
        },
        obliquity_formula: Default::default(),
    };

    let julian_day = c.julian_day();
    println!("julian day {julian_day}");

    let julian_century = c.julian_century_j2000();
    println!("julian century (J2000.0) {julian_century}");

    let mean_longitude = c.solar_geometric_mean_longitude();
    println!("mean longitude {mean_longitude}");

    let mean_anomoly = c.solar_geometric_mean_anomoly();
    println!("mean anomoly {mean_anomoly}");

    let earth_orbit_eccentricity = c.earth_orbit_eccentricity();
    println!("earth orbit eccentricity {earth_orbit_eccentricity}");

    let sun_eq_of_ctr = c.sun_equation_of_center();
    println!("sun equation of center {sun_eq_of_ctr}");

    let sun_true_long = c.sun_true_longitude();
    println!("sun true longitude {sun_true_long}");

    let sun_true_anom = c.sun_true_anomoly();
    println!("sun true anomoly {sun_true_anom}");

    let sun_rad_vec = c.sun_rad_vector_au();
    println!("sun rad vec (AU) {sun_rad_vec}");

    let sun_app_long = c.sun_apparent_longitude();
    println!("sun apparent longitude {sun_app_long}");

    let obliquity = c.obliquity();
    println!("obliquity {obliquity}");

    let sun_right_ascension = c.sun_right_ascension();
    println!("sun right ascension {sun_right_ascension}");

    let sun_declination = c.sun_declination();
    println!("sun declination {sun_declination}");

    let eot = c.equation_of_time();
    println!("equation of time {eot}");

    let ha_sunrise = c.ha_sunrise();
    println!("HA sunrise {ha_sunrise})");

    let solar_noon_utc = c.solar_noon_utc();
    println!("solar noon UTC {solar_noon_utc}");

    let sunrise = c.sunrise_time_utc();
    println!("sunrise UTC {sunrise}");

    let sunset = c.sunset_time_utc();
    println!("sunset UTC {sunset}");

    let sunlight_minutes = c.sunlight_duration_minutes();
    println!("sunlight (min) {sunlight_minutes}");

    let true_solar_time_min = c.true_solar_time_minutes();
    println!("True solar time (TT, minutes) {true_solar_time_min}");

    let hour_angle = c.hour_angle();
    println!("Hour angle {hour_angle}");

    let solar_zenith_angle = c.solar_zenith_angle();
    println!("Solar zenith angle {solar_zenith_angle}");

    let solar_elevation_angle = c.solar_elevation_angle();
    println!("Solar elevation angle {solar_elevation_angle}");

    let approx_refr = c.approx_atmospheric_refraction();
    println!("Approx atmos refraction {approx_refr}");

    let solar_elevation_corrected = c.solar_elevation_corrected();
    println!("Solar elevation corrected {solar_elevation_corrected}");

    let solar_azimuth_angle = c.solar_azimuth_angle();
    println!("Solar azimuth angle {solar_azimuth_angle}");

    let sun_vec = c.sun_vector();
    println!("Sun vector: {sun_vec:?}");
}
