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
    }
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