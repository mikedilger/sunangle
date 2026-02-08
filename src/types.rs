use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Degrees(pub f64);

impl Degrees {
    pub fn normalize(&mut self) {
        self.0 = self.0.rem_euclid(360.0);
    }

    pub fn from_radians(radians: f64) -> Self {
        Self(360.0 * radians / (2.0 * PI))
    }

    pub fn to_radians(&self) -> f64 {
        2.0 * PI * self.0 / 360.0
    }

    pub fn from_dms(degrees: i32, minutes: i32, seconds: f64) -> Degrees {
        Degrees(degrees as f64 + (minutes as f64 / 60.0) + (seconds / 60.0 / 60.0))
    }

    pub fn dms_degrees(&self) -> i8 {
        self.0 as i8
    }

    pub fn dms_minutes(&self) -> i8 {
        ((self.0).fract() * 60.0) as i8
    }

    pub fn dms_seconds(&self) -> f64 {
        ((self.0).fract() * 60.0).fract() * 60.0
    }
}

impl fmt::Display for Degrees {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}° = {}° {}’ {}”",
            self.0,
            self.dms_degrees(),
            self.dms_minutes(),
            self.dms_seconds()
        )
    }
}

impl Add for Degrees {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl Sub for Degrees {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl Mul<f64> for Degrees {
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Self(self.0 * other)
    }
}

impl Mul<Degrees> for f64 {
    type Output = Degrees;

    fn mul(self, other: Degrees) -> Self::Output {
        Degrees(self * other.0)
    }
}

impl Div<f64> for Degrees {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        Self(self.0 / other)
    }
}

pub struct Hours(pub f64);

impl Hours {
    pub fn normalize(&mut self) {
        self.0 = self.0.rem_euclid(24.0);
    }

    pub fn from_hms(hours: i32, minutes: i32, seconds: f64) -> Hours {
        Hours((hours as f64) + (minutes as f64 / 60.0) + (seconds / 60.0 / 60.0))
    }

    pub fn hours(&self) -> i32 {
        self.0 as i32
    }

    pub fn minutes(&self) -> i32 {
        (self.0.fract() * 60.0) as i32
    }

    pub fn seconds(&self) -> f64 {
        (self.0.fract() * 60.0).fract() * 60.0
    }
}

impl From<Degrees> for Hours {
    fn from(d: Degrees) -> Hours {
        Hours(d.0 / 15.0)
    }
}

impl From<Hours> for Degrees {
    fn from(h: Hours) -> Degrees {
        Degrees(h.0 * 15.0)
    }
}

impl fmt::Display for Hours {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}h = {}h {}m {}s",
            self.0,
            self.hours(),
            self.minutes(),
            self.seconds()
        )
    }
}
