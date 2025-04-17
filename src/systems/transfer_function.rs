use core::fmt;
use std::{
    any::{Any, TypeId},
    marker::PhantomData,
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use super::polynom::{Polynomial, RationalFunction};
use crate::utils::traits::{Continuous, Discrete, One, Time, Zero};

/// A transfer function representation.
///
/// This struct is parameterized by `T` (the type of coefficients) and `U` (the
/// time domain, either `Continuous` or `Discrete`).
///
/// # Type Parameters
/// - `T`: The type of the coefficients of the transfer function (e.g., `f64`,
///   `i32`, etc.).
/// - `U`: The time domain of the transfer function (either `Continuous` or
///   `Discrete`).
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
    T: Zero + One + Clone,
{
    /// Creates a new transfer function from a given numerator and denominator.
    ///
    /// # Arguments
    /// - `num`: A slice representing the numerator coefficients in ascending
    ///   order, i.e. num = num[0] + num[1]*s + ... .
    /// - `den`: A slice representing the denominator coefficients in ascending
    ///   order, i.e. den = den[0] + den[1]*s + ... .
    ///
    /// # Returns
    /// Returns a new transfer function.
    pub fn new(num: &[T], den: &[T]) -> Self {
        Self {
            rf: RationalFunction::new_from_coeffs(num, den),
            time: PhantomData::<U>,
        }
    }

    /// Creates a new transfer function from an existing rational function.
    ///
    /// # Arguments
    /// - `rf`: A `RationalFunction<T>` to be used for the transfer function.
    ///
    /// # Returns
    /// Returns a new transfer function.
    pub fn new_from_rf(rf: RationalFunction<T>) -> Self {
        Self {
            rf,
            time: PhantomData::<U>,
        }
    }
    /// Creates a new transfer function from a scalar value.
    ///
    /// # Arguments
    /// - `scalar`: A scalar value to create the transfer function.
    ///
    /// # Returns
    /// Returns a new transfer function with the scalar as both the numerator
    /// and denominator.
    pub fn new_from_scalar(scalar: T) -> Self {
        Self::new_from_rf(RationalFunction::new_from_scalar(scalar))
    }

    /// Returns the numerator polynomial of the transfer function.
    ///
    /// # Returns
    /// A slice of the numerator coefficients.
    ///
    /// # Example
    /// ```rust
    /// use control_systems_torbox::Tf;
    /// let tf = (1.0 + Tf::s())/Tf::s();
    /// let num_poly = tf.numerator();
    /// let num_coeffs = num_poly.coeffs();
    /// ```
    pub fn numerator(&self) -> &Polynomial<T> {
        self.rf.numerator()
    }

    /// Returns the denominator coefficients of the transfer function.
    ///
    /// # Returns
    /// A slice of the denominator coefficients.
    /// # Example
    /// ```rust
    /// use control_systems_torbox::Tf;
    /// let tf = (1.0 + Tf::s())/Tf::s();
    /// let den_poly = tf.denominator();
    /// let den_coeffs = den_poly.coeffs();
    /// ```
    pub fn denominator(&self) -> &Polynomial<T> {
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
    /// Normalizes the transfer function such that the highest coefficient of
    /// the denominator is one.
    ///
    /// I.e. tf = (a0 + ... + a_m s^m) / (b0 + ... + 1 s^n)
    ///
    /// # Returns
    /// A new normalized transfer function.
    pub fn normalize(&self) -> Self {
        Tf::new_from_rf(self.rf.normalize())
    }
}

impl Tf<f64, Continuous> {
    /// Returns a transfer function representing continuous-time laplace
    /// operator
    ///
    /// # Returns
    /// continuous-time lapace operator, s, `Tf<f64, Continuous>`.
    pub fn s() -> Self {
        Tf::new(&[0., 1.], &[1.])
    }
}

impl Tf<f64, Discrete> {
    /// Returns a transfer function representing discrete-time laplace operator
    ///
    /// # Returns
    /// discrete-time lapace operator, s, `Tf<f64, Discrete>`.
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
    /// Returns the degree of the numerator and denominator of the transfer
    /// function.
    ///
    /// # Returns
    /// A tuple with the degree of the numerator and the denominator.
    pub fn degree_num_den(&self) -> (usize, usize) {
        self.rf.degree_num_den()
    }

    /// Returns the relative degree of the transfer function.
    ///
    /// I.e. degree numerator - degree denominator
    ///
    /// # Returns
    /// The relative degree of the transfer function.
    pub fn relative_degree(&self) -> i32 {
        self.rf.relative_degree()
    }

    /// Returns `true` if the transfer function is proper (relative degree ≤ 0).
    ///
    /// # Returns
    /// `true` if the transfer function is proper, `false` otherwise.
    pub fn is_proper(&self) -> bool {
        self.relative_degree() <= 0
    }

    /// Returns `true` if the transfer function is strictly proper (relative
    /// degree < 0).
    ///
    /// # Returns
    /// `true` if the transfer function is strictly proper, `false` otherwise.
    pub fn is_strictly_proper(&self) -> bool {
        self.relative_degree() < 0
    }
}

macro_rules! impl_operator_tf {
    ([$(($trait:ident, $method:ident)), *]) => {
    $(
        impl<T, U> $trait for Tf<T, U>
        where
            T: $trait<Output = T> + Clone + Zero + One+ Default + Add + AddAssign + Mul<Output = T> + Neg<Output = T>,
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
    T: Neg<Output = T> + Clone + Zero + One,
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
            T: $trait<Output = T> + $assign_trait + One + Clone + Default + AddAssign + Neg<Output = T> + Mul<Output = T> + Add + Zero,
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
    /// Evaluates the transfer function at a given input `x`. Typically x = j
    /// omega
    ///
    /// # Arguments
    /// - `x`: The input value to evaluate the transfer function at.
    ///
    /// # Returns
    /// The result of evaluating the transfer function at `x`.
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
    /// Raises the transfer function to a specified integer power.
    ///
    /// # Arguments
    /// - `exp`: The exponent to which the transfer function should be raised.
    ///
    /// # Returns
    /// A new transfer function raised to the specified power.
    pub fn powi(self, exp: i32) -> Self {
        let new_rf = self.rf.powi(exp);
        Tf::new_from_rf(new_rf)
    }
}

impl<T, U: Time> Tf<T, U>
where
    T: Clone
        + Add<T, Output = T>
        + Zero
        + One
        + Default
        + Add
        + AddAssign
        + Mul<Output = T>
        + Neg<Output = T>,
{
    /// Connects two systems in **series**.
    ///
    /// Given two systems `self` and `sys2`, the output of `self` is fed as the
    /// input to `sys2`.
    ///
    /// Mathematically:
    /// ```txt
    /// u ---> [ self ] ---> [ sys2 ] ---> y
    /// ```
    ///
    /// # Returns
    /// A new system representing the series connection of `self` and `sys2`.
    pub fn series(self, sys2: Self) -> Self {
        sys2 * self
    }

    /// Connects two systems in **parallel**.
    ///
    /// The two systems `self` and `sys2` operate independently on the same
    /// input `u`, and their outputs are summed together.
    ///
    /// Mathematically:
    /// ```txt
    ///       |-----> [ self ] -----|
    /// u --->|                    (sum)---> y
    ///       |-----> [ sys2 ] -----|
    /// ```
    ///
    /// # Returns
    /// A new system representing the parallel connection of `self` and `sys2`.
    pub fn parallel(self, sys2: Self) -> Self {
        self + sys2
    }

    /// **Negative** Feedback connection of `self` with `sys2`.
    ///
    /// The system `sys2` provides feedback to `self`, forming a closed-loop
    /// system.
    ///
    /// ## Diagram:
    /// ```txt
    ///             --------        y
    /// u ---->O--->| self |--------->
    ///      -1^    --------    |
    ///        |                |
    ///        |    --------    |
    ///        -----| sys2 |<----
    ///             --------
    /// ```
    ///
    /// # Returns
    /// A new system representing the negative feedback connection of `self`
    /// with `sys2`.
    pub fn feedback(self, sys2: Self) -> Self {
        let n1 = self.numerator();
        let d1 = self.denominator();
        let n2 = sys2.numerator();
        let d2 = sys2.denominator();

        let new_num = n1.clone() * d2.clone();
        let new_den = n1.clone() * n2.clone() + d1.clone() * d2.clone();
        let rf = RationalFunction::new_from_poly(new_num, new_den);
        Tf::new_from_rf(rf)
    }
}

///////////////////////////////////////////////////////////////////////////////////
/// TESTS
///////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use crate::utils::traits::Continuous;

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

        let tf_z = Tf::<f64, Discrete>::new_from_rf(tf.rf);
        println!("Tf discrete: \n {}", tf_z);
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

        let mut tf = Tf::s();
        tf += Tf::new_from_scalar(1.0);
        let tf2 = 1.0 + Tf::s();
        assert_abs_diff_eq!(tf, tf2);
    }

    #[test]
    fn tf_math() {
        let tf = Tf::s() + 1. / (Tf::s() + 1.0);
        let tf_neg = -tf.clone();
        let tf_neg_neg = -tf_neg;
        assert_eq!(tf, tf_neg_neg);
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

    #[test]
    fn tf_interconnections() {
        let tf1 = Tf::s() / (1.0 + Tf::s());
        let tf2 = 1.0 / Tf::s();

        let tf_add = tf1.clone() + tf2.clone();
        let tf_parallel = tf1.clone().parallel(tf2.clone());
        assert_eq!(tf_add, tf_parallel);

        let tf_mul = tf2.clone() * tf1.clone();
        let tf_series = tf1.clone().series(tf2.clone());
        assert_eq!(tf_series, tf_mul);
    }
}
