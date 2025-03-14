use core::fmt;
use std::{
    any::{Any, TypeId},
    marker::PhantomData,
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use crate::{
    polynom::RationalFunction,
    traits::{Continuous, Discrete, One, Time, Zero},
};

#[derive(Clone, Debug)]
pub struct Tf<T, U: Time> {
    rf: RationalFunction<T>,
    time: PhantomData<U>,
}

impl<T, U: Time> PartialEq for Tf<T, U>
where
    T: One
        + Div<Output = T>
        + Zero
        + Clone
        + Mul<Output = T>
        + Default
        + AddAssign
        + PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.rf.eq(&other.rf)
    }
}

impl<T, U: Time> Tf<T, U>
where
    T: Zero + Clone,
{
    pub fn new(num: &[T], den: &[T]) -> Self {
        Self {
            rf: RationalFunction::new_from_coeffs(num, den),
            time: PhantomData::<U>,
        }
    }

    pub fn new_from_rf(rf: RationalFunction<T>) -> Self {
        Self {
            rf,
            time: PhantomData::<U>,
        }
    }

    pub fn numerator(&self) -> &[T] {
        self.rf.numerator()
    }

    pub fn denominator(&self) -> &[T] {
        self.rf.denominator()
    }
}

impl<T, U: Time> Tf<T, U>
where
    T: One
        + Div<Output = T>
        + Zero
        + Clone
        + Mul<Output = T>
        + Default
        + AddAssign,
{
    pub fn normalize(self) -> Self {
        Tf::new_from_rf(self.rf.normalize())
    }
}

impl Tf<f64, Continuous> {
    pub fn s() -> Self {
        Tf::new(&[0., 1.], &[1.])
    }
}

impl Tf<f64, Discrete> {
    pub fn z() -> Self {
        Tf::new(&[0., 1.], &[1.])
    }
}

impl<T, U: Time + Any> fmt::Display for Tf<T, U>
where
    T: Zero + fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let var_name = if TypeId::of::<U>() == TypeId::of::<Continuous>() {
            "s"
        } else if TypeId::of::<U>() == TypeId::of::<Discrete>() {
            "z"
        } else {
            "x"
        };
        let rf_str = self.rf.to_string_variable(var_name);
        write!(f, "{}", rf_str)
    }
}

impl<T, U: Time> Tf<T, U>
where
    T: Zero,
{
    pub fn degree_num_den(&self) -> (usize, usize) {
        self.rf.degree_num_den()
    }

    pub fn relative_degree(&self) -> i32 {
        self.rf.relative_degree()
    }

    pub fn is_proper(&self) -> bool {
        self.relative_degree() <= 0
    }

    pub fn is_strictly_proper(&self) -> bool {
        self.relative_degree() < 0
    }
}

macro_rules! impl_operator_tf {
    ([$(($trait:ident, $method:ident)), *]) => {
    $(
        impl<T, U> $trait for Tf<T, U>
        where
            T: $trait<Output = T> + Clone + Zero + Default + Add + AddAssign + Mul<Output = T> + Neg<Output = T>,
            U: Time,
        {
            type Output = Tf<T, U>;
            fn $method(self, rhs: Self) -> Self::Output {
                let new_rf = RationalFunction::$method(self.rf, rhs.rf);
                Tf::new_from_rf(new_rf)
            }
        }
    )*
    };
}
impl_operator_tf!([(Add, add), (Sub, sub), (Mul, mul), (Div, div)]);

impl<T, U> Neg for Tf<T, U>
where
    T: Neg<Output = T> + Clone + Zero,
    U: Time,
{
    type Output = Tf<T, U>;
    fn neg(self) -> Self::Output {
        let new_rf = -self.rf;
        Tf::new_from_rf(new_rf)
    }
}

macro_rules! impl_compound_assign {
    ($struct_type:ident, [$(($trait:ident, $method:ident, $assign_trait:ident, $assign_method:ident )), *]) => {
    $(
        impl<T, U: Time> $assign_trait for $struct_type<T, U>
        where
            T: $trait<Output = T> + $assign_trait + Clone + Default + AddAssign + Neg<Output = T> + Mul<Output = T> + Add + Zero,
        {
            fn $assign_method(&mut self, rhs: Self) {
                *self = self.clone().$method(rhs)
            }
        }
    )*
    };
}
impl_compound_assign!(
    Tf,
    [
        (Add, add, AddAssign, add_assign),
        (Sub, sub, SubAssign, sub_assign),
        (Mul, mul, MulAssign, mul_assign),
        (Div, div, DivAssign, div_assign)
    ]
);

macro_rules! impl_scalar_math_operator {
    ($struct_type:ident, $scalar_type:ty, [$(($operator:ident, $operator_fn:ident)), *]) => {
        $(
            impl<U: Time> $operator<$scalar_type> for $struct_type<$scalar_type, U> {
                type Output = Self;
                fn $operator_fn(self, rhs: $scalar_type) -> Self::Output {
                    let new_rf = self.rf.$operator_fn(RationalFunction::new_from_scalar(rhs));
                    Tf::new_from_rf(new_rf)
                }
            }

            impl<U: Time> $operator<$struct_type<$scalar_type, U>> for $scalar_type {
                type Output = $struct_type<$scalar_type, U>;
                fn $operator_fn(self, rhs: $struct_type<$scalar_type, U>) -> Self::Output {
                    let scalar_rf = RationalFunction::new_from_scalar(self);
                    let new_rf = scalar_rf.$operator_fn(rhs.rf);
                    Tf::new_from_rf(new_rf)
                }
            }
        )*
    };
}

impl_scalar_math_operator!(
    Tf,
    f32,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);
impl_scalar_math_operator!(
    Tf,
    f64,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);

impl_scalar_math_operator!(
    Tf,
    i8,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);
impl_scalar_math_operator!(
    Tf,
    i16,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);

impl_scalar_math_operator!(
    Tf,
    i32,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);
impl_scalar_math_operator!(
    Tf,
    i64,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);
impl_scalar_math_operator!(
    Tf,
    i128,
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);
impl_scalar_math_operator!(Tf, u8, [(Add, add), (Mul, mul), (Div, div)]);
impl_scalar_math_operator!(Tf, u16, [(Add, add), (Mul, mul), (Div, div)]);

impl_scalar_math_operator!(Tf, u32, [(Add, add), (Mul, mul), (Div, div)]);
impl_scalar_math_operator!(Tf, u64, [(Add, add), (Mul, mul), (Div, div)]);
impl_scalar_math_operator!(Tf, u128, [(Add, add), (Mul, mul), (Div, div)]);
impl<T, U: Time> Tf<T, U>
where
    T: One + Zero + Mul<Output = T> + AddAssign + Clone,
{
    pub fn eval<N>(&self, x: &N) -> N
    where
        N: Clone
            + One
            + Zero
            + Mul<N, Output = N>
            + Add<T, Output = N>
            + Div<Output = N>,
    {
        self.rf.eval(x)
    }
}

impl<T, U: Time> Tf<T, U>
where
    T: One + Clone + Zero + Add + Mul<Output = T> + AddAssign + Default,
{
    pub fn powi(self, exp: i32) -> Self {
        let new_rf = self.rf.powi(exp);
        Tf::new_from_rf(new_rf)
    }
}

///////////////////////////////////////////////////////////////////////////////////
/// TESTS
///////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use crate::traits::Continuous;

    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::Rng;

    use approx::AbsDiffEq;

    impl<U: Time> AbsDiffEq for Tf<f64, U> {
        type Epsilon = f64;

        fn default_epsilon() -> Self::Epsilon {
            1e-6
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            self.rf.abs_diff_eq(&other.rf, epsilon)
        }
    }

    #[test]
    fn eval_tf() {
        let mut rng = rand::rng();
        let tf: Tf<f64, Continuous> = Tf::new(&[1.0], &[10.0, 1.0]);
        for _ in 0..1000 {
            let x = num_complex::c64(rng.random::<f64>(), rng.random::<f64>());
            let y = tf.eval(&x);
            let y_expect = 1.0 / (10.0 + x);
            assert_abs_diff_eq!(y.re, y_expect.re);
            assert_abs_diff_eq!(y.im, y_expect.im);
        }
        println!("Tf cont: \n{}", tf);
    }

    #[test]
    fn shortcut() {
        let s = Tf::s();
        assert_eq!(s.eval(&10.), 10.);

        let z = Tf::z();
        assert_eq!(z.eval(&101.), 101.);
    }

    #[test]
    fn scalar_math() {
        let s = Tf::s();
        let sys = (2. * s + 1.) / (Tf::s() + 1.);

        let mut rng = rand::rng();
        for _ in 0..1000 {
            let x = rng.random::<f64>();
            assert_abs_diff_eq!(sys.eval(&x), (2. * x + 1.) / (x + 1.));
        }
    }

    #[test]
    fn powi_tf() {
        let sys =
            (2. * Tf::s().powi(0) + 1. * Tf::s().powi(1) + Tf::s().powi(2))
                / (1. * Tf::s().powi(1) + 2. * Tf::s().powi(2));
        let sys_new: Tf<f64, Continuous> =
            Tf::new(&[2., 1., 1.], &[0., 1., 2.]);
        assert_eq!(sys.rf, sys_new.rf);
    }

    #[test]
    fn degree() {
        let sys =
            Tf::<f64, Discrete>::new(&[0.0, 0., 0., 1., 2., 0.0], &[0.0, 1.]);
        assert_eq!(sys.degree_num_den(), (4, 1));
        assert_eq!(sys.relative_degree(), 4 - 1);
        assert_eq!(sys.is_proper(), false);
        assert_eq!(sys.is_strictly_proper(), false);

        let sys = Tf::s() / Tf::s();
        assert_eq!(sys.is_proper(), true);
        assert_eq!(sys.is_strictly_proper(), false);

        let sys = 1. / Tf::s();
        assert_eq!(sys.is_proper(), true);
        assert_eq!(sys.is_strictly_proper(), true);
    }
}
