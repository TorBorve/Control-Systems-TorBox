use std::{
    marker::PhantomData,
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use crate::{
    polynom::RationalFunction,
    traits::{One, Time, Zero},
};

#[derive(Clone, Debug)]
pub struct Tf<T, U: Time> {
    rf: RationalFunction<T>,
    time: PhantomData<U>,
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


#[cfg(test)]
mod tests {
    use crate::traits::Continuous;

    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::Rng;

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
    }
}