use std::fmt::Debug;

pub trait Zero {
    fn zero() -> Self;

    fn is_zero(&self) -> bool;
}

pub trait One {
    fn one() -> Self;

    fn is_one(&self) -> bool;
}

macro_rules! impl_traits {
    ($t:ty, $zero:expr, $one:expr) => {
        impl Zero for $t {
            fn zero() -> Self {
                $zero
            }

            fn is_zero(&self) -> bool {
                *self == $zero
            }
        }

        impl One for $t {
            fn one() -> Self {
                $one
            }

            fn is_one(&self) -> bool {
                *self == $one
            }
        }
    };
}

impl_traits!(f32, 0., 1.);
impl_traits!(f64, 0., 1.);
impl_traits!(i8, 0, 1);
impl_traits!(i16, 0, 1);
impl_traits!(i32, 0, 1);
impl_traits!(i64, 0, 1);
impl_traits!(i128, 0, 1);
impl_traits!(isize, 0, 1);
impl_traits!(u8, 0, 1);
impl_traits!(u16, 0, 1);
impl_traits!(u32, 0, 1);
impl_traits!(u64, 0, 1);
impl_traits!(u128, 0, 1);
impl_traits!(usize, 0, 1);
impl_traits!(
    num_complex::Complex32,
    num_complex::c32(0.0, 0.0),
    num_complex::c32(1.0, 0.0)
);
impl_traits!(
    num_complex::Complex64,
    num_complex::c64(0.0, 0.0),
    num_complex::c64(1.0, 0.0)
);

pub trait Mag2Db {
    fn mag2db(&self) -> Self;
}

impl Mag2Db for f64 {
    fn mag2db(&self) -> Self {
        20. * self.log10()
    }
}

impl Mag2Db for f32 {
    fn mag2db(&self) -> Self {
        20. * self.log10()
    }
}

pub trait Db2Mag {
    fn db2mag(&self) -> Self;
}

impl Db2Mag for f64 {
    fn db2mag(&self) -> Self {
        10_f64.powf(self / 20.)
    }
}

impl Db2Mag for f32 {
    fn db2mag(&self) -> Self {
        10_f32.powf(self / 20.)
    }
}

pub trait Rad2Deg {
    fn rad2deg(&self) -> Self;
}

pub trait Deg2Rad {
    fn deg2rad(&self) -> Self;
}

impl Rad2Deg for f64 {
    fn rad2deg(&self) -> Self {
        self * 180. / std::f64::consts::PI
    }
}

impl Rad2Deg for f32 {
    fn rad2deg(&self) -> Self {
        self * 180. / std::f32::consts::PI
    }
}

impl Deg2Rad for f64 {
    fn deg2rad(&self) -> Self {
        self * std::f64::consts::PI / 180.
    }
}

impl Deg2Rad for f32 {
    fn deg2rad(&self) -> Self {
        self * std::f32::consts::PI / 180.
    }
}

pub trait Time: Clone + Debug {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Continuous {}
impl Time for Continuous {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Discrete {}
impl Time for Discrete {}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI};

    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use rand::Rng;

    use super::*;

    #[test]
    fn mag_db_covert() {
        assert_abs_diff_eq!((1. as f64).mag2db(), 0.);
        assert_abs_diff_eq!((1. as f32).mag2db(), 0.);

        assert_abs_diff_eq!((20 as f64).db2mag(), 10.);
        assert_abs_diff_eq!((20 as f32).db2mag(), 10.);

        let mut rng = rand::rng();

        for _ in 0..10000 {
            let x: f64 = rng.random_range(-1.0..1000.0);
            if x < 0. {
                assert!((x as f32).mag2db().is_nan());
                assert!((x as f64).mag2db().is_nan());
            } else {
                assert_abs_diff_eq!(
                    (x as f32).mag2db().db2mag(),
                    x as f32,
                    epsilon = 1e-3
                );
                assert_abs_diff_eq!(
                    (x as f64).mag2db().db2mag(),
                    x as f64,
                    epsilon = 1e-3
                );
            }
        }
    }

    #[test]
    fn angle_convert() {
        assert_abs_diff_eq!(FRAC_PI_2.rad2deg(), 90.);
        assert_abs_diff_eq!(-180.0.deg2rad(), -PI);

        let mut rng = rand::rng();
        for _ in 0..10000 {
            let x: f64 = rng.random_range(-100.0 * PI..100.0 * PI);
            assert_relative_eq!(x.rad2deg().deg2rad(), x,);
            assert_relative_eq!((x as f32).rad2deg().deg2rad(), x as f32);
        }
    }

    #[test]
    fn one_and_zero() {
        assert!(f64::one().is_one());
        assert!(f64::zero().is_zero());
    }
}
