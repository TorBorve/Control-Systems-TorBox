use crate::traits::{One, Zero};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

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
}

impl<T> Polynomial<T>
where
    T: One + Zero + Mul<Output = T> + AddAssign + Clone,
{
    pub fn eval(&self, x: &T) -> T {
        let mut x_i = T::one();
        let mut sum = T::zero();
        for coeff_i in self.coeffs.iter() {
            sum += coeff_i.clone() * x_i.clone();
            x_i = x_i * x.clone();
        }
        sum
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
}

impl<T> RationalFunction<T>
where
    T: One + Zero + Mul<Output = T> + Div<Output = T> + AddAssign + Clone,
{
    pub fn eval(&self, x: &T) -> T {
        let num_val = self.num.eval(x);
        let den_val = self.den.eval(x);
        num_val / den_val
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::Rng;
    use std::f64;

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
    }

    #[test]
    fn rational_functions_normalize() {
        let rf =
            RationalFunction::new_from_coeffs(&[2., 0., 4., 6.], &[0., 2.0]);
        let rf = rf.normalize();
        assert_abs_diff_eq!(rf.den.coeffs.last().unwrap().clone(), 1.0);
        assert_abs_diff_eq!(rf.num.coeffs[0], 1.0);
        assert_abs_diff_eq!(rf.num.coeffs.last().unwrap().clone(), 3.0);
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
}
