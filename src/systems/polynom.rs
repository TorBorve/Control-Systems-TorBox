use crate::utils::traits::{One, Zero};
use core::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

/////////////////////////////////////////////////////////
/// Polynomial
////////////////////////////////////////////////////////
#[derive(Debug, PartialEq, Clone)]
pub struct Polynomial<T> {
    coeffs: Vec<T>,
}

impl<T: Zero + Clone> Polynomial<T> {
    pub fn new(coeffs: &[T]) -> Self {
        let poly = Self {
            coeffs: coeffs.into(),
        };
        poly.trim()
    }

    pub fn trim(mut self) -> Self {
        if let Some(pos) = self.coeffs.iter().rposition(|x| !x.is_zero()) {
            let new_len = pos + 1;
            self.coeffs.truncate(new_len);
        } else {
            self.coeffs.resize(1, T::zero());
        }

        self
    }

    pub fn new_from_scalar(value: T) -> Self {
        Self::new(&[value])
    }
}

impl<T> From<T> for Polynomial<T>
where
    T: Zero + Clone,
{
    fn from(value: T) -> Self {
        Polynomial::new_from_scalar(value)
    }
}

impl<T> Polynomial<T>
where
    T: Zero + fmt::Display,
{
    pub fn to_string_variable(&self, var_name: &str) -> String {
        let mut terms = Vec::new();
        for (idx, val) in self.coeffs.iter().rev().enumerate() {
            if val.is_zero() {
                continue;
            }

            let exp = self.coeffs.len() - idx - 1;

            let term_i = match exp {
                0 => val.to_string(),
                1 => format!("{}{}", val, var_name),
                _ => format!("{}{}^{}", val, var_name, exp),
            };
            terms.push(term_i);
        }

        terms.join(" + ").replace("+ -", "- ")
    }
}

impl<T> fmt::Display for Polynomial<T>
where
    T: Zero + fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str_poly = self.to_string_variable("x");
        write!(f, "{}", str_poly)
    }
}

impl<T> Polynomial<T>
where
    T: One + Zero + Mul<Output = T> + AddAssign + Clone,
{
    pub fn eval<N>(&self, x: &N) -> N
    where
        N: Clone + One + Zero + Mul<N, Output = N> + Add<T, Output = N>,
    {
        // a + bx + cx^2 ... = a + x(b + x(c + ...))
        self.coeffs
            .iter()
            .rev()
            .fold(N::zero(), |acc, coeff_i| acc * x.clone() + coeff_i.clone())
    }
}

impl<T> Polynomial<T>
where
    T: Zero,
{
    pub fn degree(&self) -> usize {
        let last_non_zero = self.coeffs.iter().rposition(|x| !x.is_zero());
        last_non_zero.unwrap_or(0)
    }
}

impl<T> Polynomial<T> {
    pub fn coeffs(&self) -> &[T] {
        &self.coeffs
    }
}

impl<T> Mul for Polynomial<T>
where
    T: Default + Mul<Output = T> + AddAssign + Clone,
{
    type Output = Polynomial<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result =
            vec![T::default(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for (idx_l, val_l) in self.coeffs.iter().enumerate() {
            for (idx_r, val_r) in rhs.coeffs.iter().enumerate() {
                result[idx_l + idx_r] += val_l.clone() * val_r.clone();
            }
        }

        Polynomial { coeffs: result }
    }
}

impl<T> Add for Polynomial<T>
where
    T: AddAssign + Default + Clone,
{
    type Output = Polynomial<T>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result =
            vec![T::default(); self.coeffs.len().max(rhs.coeffs.len())];
        for (idx, val) in self.coeffs.iter().enumerate() {
            result[idx] += val.clone();
        }
        for (idx, val) in rhs.coeffs.iter().enumerate() {
            result[idx] += val.clone();
        }

        Polynomial { coeffs: result }
    }
}

impl<T> Neg for Polynomial<T>
where
    T: Neg<Output = T> + Clone,
{
    type Output = Polynomial<T>;
    fn neg(self) -> Self::Output {
        Polynomial {
            coeffs: self.coeffs.into_iter().map(|x| -x).collect(),
        }
    }
}

impl<T> Sub for Polynomial<T>
where
    T: Neg<Output = T> + Clone + AddAssign + Default,
{
    type Output = Polynomial<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        let neg_rhs = -rhs;
        self + neg_rhs
    }
}

macro_rules! impl_compound_assign {
    ($struct_type:ident, [$(($trait:ident, $method:ident, $assign_trait:ident, $assign_method:ident )), *]) => {
    $(
        impl<T> $assign_trait for $struct_type<T>
        where
            T: $trait<Output = T> + $assign_trait + Clone + Default + AddAssign + Neg<Output = T> + Mul<Output = T> + Add,
        {
            fn $assign_method(&mut self, rhs: Self) {
                *self = self.clone().$method(rhs)
            }
        }
    )*
    };
}

impl_compound_assign!(
    Polynomial,
    [
        (Add, add, AddAssign, add_assign),
        (Sub, sub, SubAssign, sub_assign),
        (Mul, mul, MulAssign, mul_assign)
    ]
);

// Did this because of orphanrule. However, could also just use .into()?
macro_rules! impl_one_operator_scalar_trait {
    ($struct_type:ident, $scalar_type:ty , [$(($operator:ident, $operator_fn:ident)), *]) => {
        $(
            impl $operator<$scalar_type> for $struct_type<$scalar_type>
            {
                type Output = Self;
                fn $operator_fn(self, rhs: $scalar_type) -> Self::Output {
                    self.$operator_fn($struct_type::<$scalar_type>::new_from_scalar(rhs))
                }
            }

            impl $operator<$struct_type<$scalar_type>> for $scalar_type
            {
                type Output = $struct_type<$scalar_type>;
                fn $operator_fn(self, rhs: $struct_type<$scalar_type>) -> Self::Output {
                    let scalar = $struct_type::<$scalar_type>::new_from_scalar(self);
                    scalar.$operator_fn(rhs)
                }
            }
        )*
    };
}

impl_one_operator_scalar_trait!(
    Polynomial,
    f64,
    [(Sub, sub), (Add, add), (Mul, mul)]
);

impl_one_operator_scalar_trait!(
    Polynomial,
    f32,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(
    Polynomial,
    i8,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(
    Polynomial,
    i16,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(
    Polynomial,
    i32,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(
    Polynomial,
    i64,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(
    Polynomial,
    i128,
    [(Sub, sub), (Add, add), (Mul, mul)]
);

impl_one_operator_scalar_trait!(
    Polynomial,
    isize,
    [(Sub, sub), (Add, add), (Mul, mul)]
);
impl_one_operator_scalar_trait!(Polynomial, u8, [(Add, add), (Mul, mul)]);
impl_one_operator_scalar_trait!(Polynomial, u16, [(Add, add), (Mul, mul)]);

impl_one_operator_scalar_trait!(Polynomial, u32, [(Add, add), (Mul, mul)]);
impl_one_operator_scalar_trait!(Polynomial, u64, [(Add, add), (Mul, mul)]);
impl_one_operator_scalar_trait!(Polynomial, u128, [(Add, add), (Mul, mul)]);
impl_one_operator_scalar_trait!(Polynomial, usize, [(Add, add), (Mul, mul)]);

/////////////////////////////////////////////////////////////
/// Rational function
////////////////////////////////////////////////////////////
#[derive(Debug, Clone)]
pub struct RationalFunction<T> {
    num: Polynomial<T>,
    den: Polynomial<T>,
}

impl<T: Zero + Clone> RationalFunction<T> {
    pub fn new_from_coeffs(num: &[T], den: &[T]) -> Self {
        let num = Polynomial::new(num);
        let den = Polynomial::new(den);
        Self { num, den }
    }

    pub fn new_from_poly(
        num_poly: Polynomial<T>,
        den_poly: Polynomial<T>,
    ) -> Self {
        Self {
            num: num_poly,
            den: den_poly,
        }
    }
}

impl<T> RationalFunction<T> {
    pub fn numerator(&self) -> &Polynomial<T> {
        &self.num
    }

    pub fn denominator(&self) -> &Polynomial<T> {
        &self.den
    }
}

impl<T> RationalFunction<T>
where
    T: Zero + Clone + One,
{
    pub fn new_from_scalar(scalar: T) -> Self {
        Self::new_from_coeffs(&[scalar], &[T::one()])
    }
}

impl<T> From<T> for RationalFunction<T>
where
    T: Zero + Clone + One,
{
    fn from(value: T) -> Self {
        Self::new_from_scalar(value)
    }
}

impl<T> RationalFunction<T>
where
    T: Zero + fmt::Display,
{
    pub fn to_string_variable(&self, var_name: &str) -> String {
        let num_str = self.num.to_string_variable(var_name);
        let den_str = self.den.to_string_variable(var_name);

        let dash_line = "\u{2500}".repeat(num_str.len().max(den_str.len()));

        num_str + "\n" + &dash_line + "\n" + &den_str
    }
}

impl<T> fmt::Display for RationalFunction<T>
where
    T: Zero + fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string_variable("x"))
    }
}

impl<T> RationalFunction<T>
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
        let num_val = self.num.eval(x);
        let den_val = self.den.eval(x);
        num_val / den_val
    }
}

impl<T> RationalFunction<T>
where
    T: Zero,
{
    pub fn degree_num_den(&self) -> (usize, usize) {
        (self.num.degree(), self.den.degree())
    }

    pub fn relative_degree(&self) -> i32 {
        self.num.degree() as i32 - self.den.degree() as i32
    }
}

impl<T> RationalFunction<T>
where
    T: One
        + Div<Output = T>
        + Zero
        + Clone
        + Mul<Output = T>
        + Default
        + AddAssign,
{
    pub fn normalize(&self) -> Self {
        let highest_order_gain =
            self.den.coeffs.last().unwrap_or(&T::one()).clone();

        let inv_highest_order_gain =
            Polynomial::new(&[T::one() / highest_order_gain]);
        let new_den = inv_highest_order_gain.clone() * self.den.clone();
        let new_num = inv_highest_order_gain * self.num.clone();

        Self {
            num: new_num,
            den: new_den,
        }
    }
}

impl<T> RationalFunction<T>
where
    T: One + Clone + Zero + Add + Mul<Output = T> + AddAssign + Default,
{
    pub fn powi(self, exp: i32) -> Self {
        let base = self;
        let mut result = RationalFunction::new_from_scalar(T::one());
        for _ in 0..exp.abs() {
            result = result * base.clone();
        }
        if exp < 0 {
            result = RationalFunction::new_from_scalar(T::one()) / result
        }
        result
    }
}

impl<T> PartialEq for RationalFunction<T>
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
        let lhs = self.normalize();
        let rhs = other.normalize();
        lhs.num == rhs.num && lhs.den == rhs.den
    }
}

impl<T> Mul for RationalFunction<T>
where
    T: Add + Default + Mul<Output = T> + AddAssign + Clone,
{
    type Output = RationalFunction<T>;
    fn mul(self, rhs: Self) -> Self::Output {
        let new_num = self.num * rhs.num;
        let new_den = self.den * rhs.den;

        Self {
            num: new_num,
            den: new_den,
        }
    }
}

impl<T> Div for RationalFunction<T>
where
    T: Clone + Default + Mul<Output = T> + AddAssign,
{
    type Output = RationalFunction<T>;
    fn div(self, rhs: Self) -> Self::Output {
        let new_num = self.num * rhs.den;
        let new_den = self.den * rhs.num;

        Self {
            num: new_num,
            den: new_den,
        }
    }
}

impl<T> Add for RationalFunction<T>
where
    T: AddAssign + Default + Clone + Mul<Output = T>,
{
    type Output = RationalFunction<T>;
    fn add(self, rhs: Self) -> Self::Output {
        let new_den = self.den.clone() * rhs.den.clone();
        let new_num = self.num * rhs.den + rhs.num * self.den;

        Self {
            num: new_num,
            den: new_den,
        }
    }
}

impl<T> Neg for RationalFunction<T>
where
    T: Neg<Output = T> + Clone,
{
    type Output = RationalFunction<T>;
    fn neg(self) -> Self::Output {
        let new_num = -self.num;

        Self {
            num: new_num,
            den: self.den,
        }
    }
}

impl<T> Sub for RationalFunction<T>
where
    T: Neg<Output = T> + Clone + AddAssign + Default + Mul<Output = T>,
{
    type Output = RationalFunction<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl_compound_assign!(
    RationalFunction,
    [
        (Add, add, AddAssign, add_assign),
        (Sub, sub, SubAssign, sub_assign),
        (Mul, mul, MulAssign, mul_assign),
        (Div, div, DivAssign, div_assign)
    ]
);

impl_one_operator_scalar_trait!(
    RationalFunction,
    f64,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);

impl_one_operator_scalar_trait!(
    RationalFunction,
    f32,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    i8,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    i16,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    i32,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    i64,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    i128,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);

impl_one_operator_scalar_trait!(
    RationalFunction,
    isize,
    [(Sub, sub), (Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    u8,
    [(Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    u16,
    [(Add, add), (Mul, mul), (Div, div)]
);

impl_one_operator_scalar_trait!(
    RationalFunction,
    u32,
    [(Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    u64,
    [(Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    u128,
    [(Add, add), (Mul, mul), (Div, div)]
);
impl_one_operator_scalar_trait!(
    RationalFunction,
    usize,
    [(Add, add), (Mul, mul), (Div, div)]
);

/////////////////////////////////////////////////////////////////////////////////////////////
/// Tests
//////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::Rng;
    use std::f64;

    use approx::AbsDiffEq;

    impl AbsDiffEq for Polynomial<f64> {
        type Epsilon = f64;

        fn default_epsilon() -> Self::Epsilon {
            1e-6
        }
        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            let max_len = self.coeffs.len().max(other.coeffs.len());

            let mut lhs_vals = self.coeffs.to_vec();
            lhs_vals.resize(max_len, 0.);
            let mut rhs_vals = other.coeffs.to_vec();
            rhs_vals.resize(max_len, 0.);

            lhs_vals
                .iter()
                .zip(rhs_vals.iter())
                .all(|(x, y)| f64::abs_diff_eq(&x, &y, epsilon))
        }
    }

    impl AbsDiffEq for RationalFunction<f64> {
        type Epsilon = f64;

        fn default_epsilon() -> Self::Epsilon {
            1e-6
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            let lhs = self.normalize();
            let rhs = other.normalize();

            lhs.num.abs_diff_eq(&rhs.num, epsilon)
                && lhs.den.abs_diff_eq(&rhs.den, epsilon)
        }
    }

    #[test]
    fn polynomial() {
        let mut rng = rand::rng();

        let mut gen_rand_coeffs = || {
            let rand_coeffs: Vec<f64> =
                (0..10).map(|_| rng.random_range(1..=100) as f64).collect();
            rand_coeffs
        };

        let p = Polynomial::new(&[1., 2., 0.]);
        assert_abs_diff_eq!(p.eval(&1.), 1. + 2. * 1., epsilon = 1e-3);
        assert_eq!(p.coeffs.len(), 2); // Trim zero

        {
            let p1 = Polynomial::new(gen_rand_coeffs().as_slice());
            let x = -1.90;
            let p1_val = p1.eval(&x);

            let p2 = Polynomial::new(gen_rand_coeffs().as_slice());
            let p2_val = p2.eval(&x);

            let p_mult = p1.clone() * p2.clone();
            assert_abs_diff_eq!(
                p_mult.eval(&x),
                p2_val * p1_val,
                epsilon = 1e-3
            );
            let p_add = p1.clone() + p2.clone();
            assert_abs_diff_eq!(
                p_add.eval(&x),
                p2_val + p1_val,
                epsilon = 1e-3
            );
            let p_sub = p1.clone() - p2.clone();
            assert_abs_diff_eq!(
                p_sub.eval(&x),
                p1_val - p2_val,
                epsilon = 1e-3
            );
        }

        {
            println!("poly: {}", p);
            let p_str = p.to_string();
            println!("poly to string: {}", p_str);
        }
    }

    #[test]
    fn polynomial_assignment() {
        let mut rng = rand::rng();
        let mut gen_rand_coeffs = || {
            let rand_coeffs: Vec<f64> =
                (0..10).map(|_| rng.random_range(1..=100) as f64).collect();
            rand_coeffs
        };
        let p1 = Polynomial::new(gen_rand_coeffs().as_slice());
        let p2 = Polynomial::new(gen_rand_coeffs().as_slice());

        let mut p = p1.clone();
        p += p2.clone();
        assert_eq!(p, p1.clone() + p2.clone());

        let mut p = p1.clone();
        p -= p2.clone();
        assert_eq!(p, p1.clone() - p2.clone());

        let mut p = p1.clone();
        p *= p2.clone();
        assert_eq!(p, p1 * p2);
    }

    #[test]
    fn polynomial_scalar() {
        let p = Polynomial::new_from_scalar(1.0);
        assert_eq!(p.clone() - 1.0, 0.0.into());
        assert_eq!(p.clone() + 1.0, 2.0.into());
        assert_eq!(p.clone() * 2.0, 2.0.into());
        assert_eq!(p.clone() + 10., 10. + p.clone());
    }

    #[test]
    fn rational_functions_normalize() {
        let rf =
            RationalFunction::new_from_coeffs(&[2., 0., 4., 6.], &[0., 2.0]);
        let rf = rf.normalize();
        assert_abs_diff_eq!(rf.den.coeffs.last().unwrap().clone(), 1.0);
        assert_abs_diff_eq!(rf.num.coeffs[0], 1.0);
        assert_abs_diff_eq!(rf.num.coeffs.last().unwrap().clone(), 3.0);

        println!("Rational function normalize: \n{}", rf);
    }

    #[test]
    fn rational_functions_eval() {
        let mut rng = rand::rng();
        let coeffs: Vec<f64> = (0..100).map(|_| rng.random::<f64>()).collect();
        let rf = RationalFunction::new_from_coeffs(
            coeffs.as_slice(),
            coeffs.as_slice(),
        );

        for x in (0..1000).map(|_| rng.random::<f64>()) {
            assert_abs_diff_eq!(rf.eval(&x), 1.0);
        }
    }

    #[test]
    fn rational_functions_math() {
        let mut rng = rand::rng();
        for _ in 0..10 {
            let num: Vec<f64> = (0..10).map(|_| rng.random::<f64>()).collect();
            let den: Vec<f64> = (0..10).map(|_| rng.random::<f64>()).collect();
            let rf = RationalFunction::new_from_coeffs(
                num.as_slice(),
                den.as_slice(),
            );

            for x in (0..100).map(|_| rng.random::<f64>()) {
                let zero = rf.clone() - rf.clone();
                assert_abs_diff_eq!(zero.eval(&x), 0.0, epsilon = 1e-6);
                let unit = (rf.clone() * rf.clone()) / rf.clone() / rf.clone();
                assert_abs_diff_eq!(unit.eval(&x), 1.0, epsilon = 1e-6);
                let zero = rf.clone() + rf.clone() - rf.clone() - rf.clone();
                assert_abs_diff_eq!(zero.eval(&x), 0.0, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn rational_functions_assignment() {
        let rf1 =
            RationalFunction::new_from_coeffs(&[16.0, 32.0, 1.0], &[2.0, 8.0]);
        let rf2 =
            RationalFunction::new_from_coeffs(&[1.0, 2.0], &[8.0, 4.0, 2.0]);
        let mut rf = rf1.clone();
        rf += rf2.clone();
        assert_eq!(rf, rf1.clone() + rf2.clone());

        let mut rf = rf1.clone();
        rf -= rf2.clone();
        assert_eq!(rf, rf1.clone() - rf2.clone());

        let mut rf = rf1.clone();
        rf *= rf2.clone();
        assert_eq!(rf, rf1.clone() * rf2.clone());

        let mut rf = rf1.clone();
        rf /= rf2.clone();
        assert_eq!(rf, rf1.clone() / rf2.clone());
    }

    #[test]
    fn rational_functions_scalar_math() {
        let rf: RationalFunction<f64> = 1.0.into();
        assert_eq!(rf.clone() - 1.0, 0.0.into());
        assert_eq!(rf.clone() + 1.0, 2.0.into());
        assert_eq!(rf.clone() * 2.0, 2.0.into());
        assert_eq!(rf.clone() / 2.0, 0.5.into());
    }

    #[test]
    fn rational_functions_powi() {
        let rf: RationalFunction<f64> = 1.0.into();
        for i in 0..100 {
            let new_rf = rf.clone().powi(i);
            assert_eq!(new_rf, rf);
        }

        let rf =
            RationalFunction::new_from_coeffs(&[1., 2., 4.], &[1., 4., 2., 8.]);
        assert_eq!(rf.clone().powi(0), RationalFunction::new_from_scalar(1.));
        assert_eq!(rf.clone().powi(1), rf.clone());
        assert_eq!(rf.clone().powi(2), rf.clone() * rf.clone());
        assert_eq!(rf.clone().powi(3), rf.clone() * rf.clone() * rf.clone());
    }

    #[test]
    fn degree() {
        let pol = Polynomial::new(&[0.0, 1.0, 2.0, 0.0, 0.0]);
        assert_eq!(pol.degree(), 2);

        let pol2 = pol.clone();
        let rf = RationalFunction::new_from_poly(pol, pol2);

        assert_eq!(rf.relative_degree(), 0);
        assert_eq!(rf.degree_num_den(), (2, 2));
    }

    #[test]
    fn abs_diff_eq_test() {
        assert!(Polynomial::<f64>::default_epsilon() > 0.);
        assert!(RationalFunction::<f64>::default_epsilon() > 0.);

        let pol = Polynomial::new(&[0.0, 1.0, 2.0, 0.0, 0.0]);
        assert_abs_diff_eq!(pol.clone(), pol);

        let rf =
            RationalFunction::new_from_coeffs(&[1., 2., 4.], &[1., 4., 2., 8.]);
        assert_abs_diff_eq!(rf.clone(), rf);
    }
}
