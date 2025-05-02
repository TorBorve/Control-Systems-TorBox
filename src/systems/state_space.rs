use std::{
    any::TypeId,
    error::Error,
    fmt,
    marker::PhantomData,
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use crate::{Continuous, Discrete, utils::traits::Time};
extern crate nalgebra as na;
use na::DMatrix;

/// A state-space representation of a system.
///
/// # Members
///
/// - `a`: State Matrix
/// - `b`: Input Matrix
/// - `c`: Output Matrix
/// - `d`: Feedthrough Matrix
///
/// # Type Parameters
/// - `U`: The time domain. `Continuous` or `Discrete`
#[derive(Debug, Clone, PartialEq)]
pub struct Ss<U: Time> {
    a: DMatrix<f64>,
    b: DMatrix<f64>,
    c: DMatrix<f64>,
    d: DMatrix<f64>,
    time: PhantomData<U>,
}

impl<U: Time + 'static> Ss<U> {
    /// Creates a new state-space system from the given matrices.
    ///
    /// # Parameters:
    /// - `a`: The state matrix (A) must be square.
    /// - `b`: The input matrix (B) must have rows equal to the number of
    ///   states.
    /// - `c`: The output matrix (C) must have columns equal to the number of
    ///   states.
    /// - `d`: The feedthrough matrix (D) must have rows equal to the number of
    ///   outputs and columns equal to the number of inputs.
    ///
    /// # Returns:
    /// - `Ok(Self)`: A state-space system on success.
    /// - `Err`: An error if the matrix dimensions do not match the expected
    ///   sizes.
    ///
    /// # Errors:
    /// Returns an error if the dimensions of the matrices are inconsistent.
    pub fn new(
        a: DMatrix<f64>,
        b: DMatrix<f64>,
        c: DMatrix<f64>,
        d: DMatrix<f64>,
    ) -> Result<Self, Box<dyn Error + 'static>> {
        let ss = Self {
            a,
            b,
            c,
            d,
            time: PhantomData::<U>,
        };
        ss.is_valid()?;
        Ok(ss)
    }

    /// Returns a reference to the state matrix (A).
    pub fn a(&self) -> &DMatrix<f64> {
        self.assert_valid();
        &self.a
    }

    /// Returns a reference to the state matrix (B).
    pub fn b(&self) -> &DMatrix<f64> {
        self.assert_valid();
        &self.b
    }

    /// Returns a reference to the state matrix (C).
    pub fn c(&self) -> &DMatrix<f64> {
        self.assert_valid();
        &self.c
    }

    /// Returns a reference to the state matrix (D).
    pub fn d(&self) -> &DMatrix<f64> {
        self.assert_valid();
        &self.d
    }

    /// Returns the number of outputs in the system.
    ///
    /// # Returns:
    /// - The number of rows in the output matrix (C), which equals the number
    ///   of outputs.
    ///
    /// # Panics:
    /// Panics if the number of rows in C does not match the number of rows in
    /// D.
    pub fn noutputs(&self) -> usize {
        self.assert_valid();
        self.c.nrows()
    }

    /// Returns the number of inputs in the system.
    ///
    /// # Returns:
    /// - The number of columns in the input matrix (B), which equals the number
    ///   of inputs.
    ///
    /// # Panics:
    /// Panics if the number of columns in D does not match the number of
    /// columns in B.
    pub fn ninputs(&self) -> usize {
        self.assert_valid();
        self.b.ncols()
    }

    /// Returns the order of the system, which is the number of states.
    ///
    /// # Returns:
    /// - The number of rows (or columns) in the state matrix (A), representing
    ///   the system order.
    ///
    /// # Panics:
    /// Panics if the state matrix (A) is not square.
    pub fn order(&self) -> usize {
        self.assert_valid();
        self.a.nrows()
    }

    /// Returns the shape of the system as a tuple (number of outputs, number of
    /// inputs).
    ///
    /// # Returns:
    /// - A tuple `(ny, nu)` where `ny` is the number of outputs and `nu` is the
    ///   number of inputs.
    pub fn shape(&self) -> (usize, usize) {
        self.assert_valid();
        (self.noutputs(), self.ninputs())
    }

    /// Verifies the consistency of the matrix dimensions for state-space
    /// representation.
    ///
    /// # Parameters:
    /// - `a`: The state matrix (A), must be square.
    /// - `b`: The input matrix (B), must have the same number of rows as A.
    /// - `c`: The output matrix (C), must have the same number of columns as A.
    /// - `d`: The feedthrough matrix (D), must have the same number of rows as
    ///   C and columns as B.
    ///
    /// # Returns:
    /// - `Ok(())`: If all matrix dimensions are consistent.
    /// - `Err`: If any matrix dimensions are inconsistent, an error is returned
    ///   with a descriptive message.
    pub fn verify_dimensions(
        a: &DMatrix<f64>,
        b: &DMatrix<f64>,
        c: &DMatrix<f64>,
        d: &DMatrix<f64>,
    ) -> Result<(), std::io::Error> {
        if !a.is_square() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "A matrix must be square",
            ));
        }
        let nx = a.nrows();
        if b.nrows() != nx {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Dimensions of A and B is not consisten",
            ));
        }
        let nu = b.ncols();
        if c.ncols() != nx {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Dimensions of A and C is not consistent",
            ));
        }
        let ny = c.nrows();

        if d.ncols() != nu || d.nrows() != ny {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Dimensions of B, C, D is not consistent",
            ));
        }
        Ok(())
    }

    fn is_valid(&self) -> Result<(), String> {
        Ss::<U>::verify_dimensions(&self.a, &self.b, &self.c, &self.d)
            .map_err(|e| e.to_string())
    }

    fn assert_valid(&self) {
        let is_valid = self.is_valid();
        assert!(
            is_valid.is_ok(),
            "Ss is invalid: Message: {}, Ss: {}",
            is_valid.unwrap_err(),
            self
        );
    }

    /// Creates a new state-space system representing a static gain.
    ///
    /// This function constructs a state-space system that consists only of a
    /// gain matrix, without any internal states.
    ///
    /// # Parameters:
    /// - `gain`: The scalar gain value.
    ///
    /// # Returns:
    /// - A state-space system representing `y = gain * u`, where `y` is the
    ///   output and `u` is the input.
    ///
    /// # Panics:
    /// Panics if the system creation fails, which should not happen for valid
    /// numerical inputs.
    pub fn new_from_scalar(gain: f64) -> Self {
        Self::new(
            DMatrix::zeros(0, 0),
            DMatrix::zeros(0, 1),
            DMatrix::zeros(1, 0),
            gain * DMatrix::identity(1, 1),
        )
        .unwrap()
    }

    /// Computes the inverse of a state-space system, if possible.
    ///
    /// # Warning
    /// The need for inverse often arrised due to feedback connections. Such as
    /// T = L/(1 + L). For such cases use the `feedback` function, which is more
    /// numerically stable.
    ///
    /// # Errors
    /// Returns an error if the direct transmission matrix (D) is not
    /// invertible.
    pub fn inv(&self) -> Result<Self, Box<dyn Error + 'static>> {
        assert_eq!(
            self.ninputs(),
            self.noutputs(),
            "System must be square to find the inverse"
        );
        self.assert_valid();
        let d_inv = self
            .d()
            .clone()
            .try_inverse()
            .ok_or("Matrix D is not invertible")?;

        let a_new = self.a() - self.b() * &d_inv * self.c();
        let b_new = self.b() * &d_inv;
        let c_new = -&d_inv * self.c();
        let d_new = d_inv;

        Ss::<U>::new(a_new, b_new, c_new, d_new)
    }

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
    pub fn series(&self, sys2: &Self) -> Self {
        assert_eq!(
            self.noutputs(),
            sys2.ninputs(),
            "Self output must be same as sys2 inputs"
        );
        self.assert_valid();
        sys2.assert_valid();

        let nx2 = sys2.order();
        let nu2 = sys2.ninputs();
        let ny2 = sys2.noutputs();
        let nx1 = self.order();
        let nu1 = self.ninputs();
        let ny1 = self.noutputs();
        assert_eq!(nu2, ny1);

        let mut a_new = DMatrix::zeros(nx1 + nx2, nx1 + nx2);
        a_new.view_mut((0, 0), (nx1, nx1)).copy_from(self.a());
        a_new
            .view_mut((nx1, 0), (nx2, nx1))
            .copy_from(&(sys2.b() * self.c()));
        a_new.view_mut((nx1, nx1), (nx2, nx2)).copy_from(sys2.a());

        let mut b_new = DMatrix::zeros(nx1 + nx2, nu1);
        b_new.view_mut((0, 0), (nx1, nu1)).copy_from(self.b());
        b_new
            .view_mut((nx1, 0), (nx2, nu1))
            .copy_from(&(sys2.b() * self.d()));

        let mut c_new = DMatrix::zeros(ny2, nx1 + nx2);
        c_new
            .view_mut((0, 0), (ny2, nx1))
            .copy_from(&(sys2.d() * self.c()));
        c_new.view_mut((0, nx1), (ny2, nx2)).copy_from(sys2.c());

        let d_new = sys2.d() * self.d();

        Ss::<U>::new(a_new, b_new, c_new, d_new).unwrap()
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
    pub fn parallel(&self, sys2: &Self) -> Self {
        assert_eq!(
            self.shape(),
            sys2.shape(),
            "System must have same shape to connnect in parallel"
        );
        self.assert_valid();
        sys2.assert_valid();

        let nu = self.ninputs();
        let ny = self.noutputs();

        let nx1 = self.order();
        let nx2 = sys2.order();

        let mut a_new = DMatrix::zeros(nx1 + nx2, nx1 + nx2);
        a_new.view_mut((0, 0), (nx1, nx1)).copy_from(self.a());
        a_new.view_mut((nx1, nx1), (nx2, nx2)).copy_from(sys2.a());

        let mut b_new = DMatrix::zeros(nx1 + nx2, nu);
        b_new.view_mut((0, 0), (nx1, nu)).copy_from(self.b());
        b_new.view_mut((nx1, 0), (nx2, nu)).copy_from(sys2.b());

        let mut c_new = DMatrix::zeros(ny, nx1 + nx2);
        c_new.view_mut((0, 0), (ny, nx1)).copy_from(self.c());

        c_new.view_mut((0, nx1), (ny, nx2)).copy_from(sys2.c());
        let d_new = self.d() + sys2.d();

        Ss::<U>::new(a_new, b_new, c_new, d_new).unwrap()
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
    /// # Assumptions
    /// - The number of inputs of `self` must match the number of outputs of
    ///   `sys2`.
    /// - The number of outputs of `self` must match the number of inputs of
    ///   `sys2`.
    ///
    /// # Returns
    /// A new system representing the negative feedback connection of `self`
    /// with `sys2`.
    ///
    /// # Panics
    /// This function will panic if the input-output dimensions of `self` and
    /// `sys2` do not match.
    pub fn feedback(self, sys2: Self) -> Self {
        assert_eq!(self.ninputs(), sys2.noutputs());
        assert_eq!(self.noutputs(), sys2.ninputs());
        self.assert_valid();
        sys2.assert_valid();

        let i_d1d2 = DMatrix::identity(self.noutputs(), self.noutputs())
            + self.d() * sys2.d();

        let i_d1d2_inv = i_d1d2.try_inverse().unwrap();

        let mut c_new =
            DMatrix::zeros(self.noutputs(), self.order() + sys2.order());
        c_new
            .view_mut((0, 0), (self.noutputs(), self.order()))
            .copy_from(self.c());
        c_new
            .view_mut((0, self.order()), (self.noutputs(), sys2.order()))
            .copy_from(&(-self.d() * sys2.c()));
        c_new = i_d1d2_inv.clone() * c_new;

        let d_new = i_d1d2_inv * self.d();

        let mut a1_new =
            DMatrix::zeros(self.order(), self.order() + sys2.order());
        a1_new
            .view_mut((0, 0), (self.order(), self.order()))
            .copy_from(self.a());
        a1_new
            .view_mut((0, self.order()), (self.order(), sys2.order()))
            .copy_from(&(-self.b() * sys2.c()));
        a1_new -= self.b() * sys2.d() * &c_new;

        let mut a2_new =
            DMatrix::zeros(sys2.order(), self.order() + sys2.order());
        a2_new
            .view_mut((0, self.order()), (sys2.order(), sys2.order()))
            .copy_from(sys2.a());
        a2_new += sys2.b() * &c_new;

        let mut a_new = DMatrix::zeros(
            self.order() + sys2.order(),
            self.order() + sys2.order(),
        );
        a_new
            .view_mut((0, 0), (self.order(), self.order() + sys2.order()))
            .copy_from(&a1_new);
        a_new
            .view_mut(
                (self.order(), 0),
                (sys2.order(), self.order() + sys2.order()),
            )
            .copy_from(&a2_new);

        let mut b_new =
            DMatrix::zeros(self.order() + sys2.order(), self.ninputs());
        b_new
            .view_mut((0, 0), (self.order(), self.ninputs()))
            .copy_from(&(self.b() - self.b() * sys2.d() * &d_new));
        b_new
            .view_mut((self.order(), 0), (sys2.order(), self.ninputs()))
            .copy_from(&(sys2.b() * &d_new));

        Ss::<U>::new(a_new, b_new, c_new, d_new).unwrap()
    }
}

impl<'b, U: Time + 'static> Add<&'b Ss<U>> for &Ss<U> {
    type Output = Ss<U>;

    /// Adds two state-space systems in parallel.
    ///
    /// y = y1 + y2
    ///
    /// # Panics
    /// Panics if the input-output dimensions of both systems do not match.
    fn add(self, rhs: &'b Ss<U>) -> Self::Output {
        assert_eq!(self.shape(), rhs.shape());
        self.assert_valid();
        rhs.assert_valid();

        self.parallel(rhs)
    }
}

impl<U: Time + 'static> Add for Ss<U> {
    type Output = Ss<U>;

    /// Adds two state-space systems in parallel.
    ///
    /// y = y1 + y2
    ///
    /// # Panics
    /// Panics if the input-output dimensions of both systems do not match.
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<U: Time + 'static> Neg for Ss<U> {
    type Output = Ss<U>;

    /// Negates the output of a state-space system.
    ///
    /// I.e. y_new = -y    
    fn neg(mut self) -> Self::Output {
        self.c *= -1.0;
        self.d *= -1.0;
        self.assert_valid();
        self
    }
}

impl<'b, U: Time + 'static> Sub<&'b Ss<U>> for &Ss<U> {
    type Output = Ss<U>;

    /// Subtracts two state-space systems by adding one and negating the other.
    fn sub(self, rhs: &'b Ss<U>) -> Self::Output {
        self.assert_valid();
        rhs.assert_valid();

        // TODO: should not need to use clone
        self + &(-rhs.clone())
    }
}

impl<U: Time + 'static> Sub for Ss<U> {
    type Output = Ss<U>;

    /// Subtracts two state-space systems by adding one and negating the other.
    fn sub(self, rhs: Self) -> Self::Output {
        self.assert_valid();
        rhs.assert_valid();
        self + (-rhs)
    }
}

impl<'b, U: Time + 'static> Mul<&'b Ss<U>> for &Ss<U> {
    type Output = Ss<U>;

    /// Multiplies two state-space systems in series (cascading connection).
    ///
    /// The resulting system represents `y <- self <- rhs <- u`, where `self`
    /// follows `rhs`.
    ///
    /// # Panics
    /// Panics if the output size of `rhs` does not match the input size of
    /// `self`.
    fn mul(self, rhs: &'b Ss<U>) -> Self::Output {
        self.assert_valid();
        rhs.assert_valid();

        rhs.series(self)
    }
}

impl<U: Time + 'static> Mul for Ss<U> {
    type Output = Ss<U>;

    /// Multiplies two state-space systems in series (cascading connection).
    ///
    /// The resulting system represents `y <- self <- rhs <- u`, where `self`
    /// follows `rhs`.
    ///
    /// # Panics
    /// Panics if the output size of `rhs` does not match the input size of
    /// `self`.
    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<'b, U: Time + 'static> Div<&'b Ss<U>> for &Ss<U> {
    type Output = Ss<U>;

    /// Divides one state-space system by another (right multiplication by
    /// inverse).
    ///
    /// # Warning
    /// The need for division often arrised due to feedback connections. Such as
    /// T = L/(1 + L). For such cases use the `feedback` function, which is more
    /// numerically stable.
    ///
    /// # Panics
    /// Panics if the inverse of `rhs` cannot be computed.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &'b Ss<U>) -> Self::Output {
        self.assert_valid();
        rhs.assert_valid();

        self * &rhs.inv().unwrap()
    }
}

impl<U: Time + 'static> Div for Ss<U> {
    type Output = Ss<U>;

    /// Divides one state-space system by another (right multiplication by
    /// inverse).
    ///
    /// # Warning
    /// The need for division often arrised due to feedback connections. Such as
    /// T = L/(1 + L). For such cases use the `feedback` function, which is more
    /// numerically stable.
    ///
    /// # Panics
    /// Panics if the inverse of `rhs` cannot be computed.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}

/// Implements compound assignment operators (`+=`, `-=`, `*=`, `/=`) for
/// state-space systems.
macro_rules! impl_compound_assign {
    ($struct_type:ident, [$(($trait:ident, $method:ident, $assign_trait:ident, $assign_method:ident )), *]) => {
    $(
        impl<U: Time + 'static> $assign_trait for $struct_type<U>
        {
            fn $assign_method(&mut self, rhs: Self) {
                *self = self.clone().$method(rhs)
            }
        }

        impl<'a, U: Time + 'static> $assign_trait<&'a $struct_type<U>> for $struct_type<U>
        {
            fn $assign_method(&mut self, rhs: &'a $struct_type<U>) {
                *self = <&Self as $trait<&Self>>::$method(self, rhs);
            }
        }

    )*
    };
}
impl_compound_assign!(
    Ss,
    [
        (Add, add, AddAssign, add_assign),
        (Sub, sub, SubAssign, sub_assign),
        (Mul, mul, MulAssign, mul_assign),
        (Div, div, DivAssign, div_assign)
    ]
);

/// Implements scalar multiplication and division for state-space systems.
macro_rules! impl_scalar_math_operator_ss {
    ([$(($operator:ident, $operator_fn:ident)), *]) => {
        $(
            impl<U: Time + 'static> $operator<f64> for Ss<U> {
                type Output = Self;
                fn $operator_fn(self, rhs: f64) -> Self::Output {
                    let rhs_ss = Ss::<U>::new_from_scalar(rhs);
                    self.$operator_fn(rhs_ss)
                }
            }

            impl<U: Time + 'static> $operator<Ss<U>> for f64 {
                type Output = Ss<U>;
                fn $operator_fn(self, rhs: Ss<U>) -> Self::Output {
                    let lhs_ss = Ss::<U>::new_from_scalar(self);
                    lhs_ss.$operator_fn(rhs)
                }
            }

            impl<U: Time + 'static> $operator<f64> for &Ss<U> {
                type Output = Ss<U>;
                fn $operator_fn(self, rhs: f64) -> Self::Output {
                    let rhs_ss = Ss::<U>::new_from_scalar(rhs);
                    self.$operator_fn(&rhs_ss)
                }
            }

            impl<U: Time + 'static> $operator<&Ss<U>> for f64 {
                type Output = Ss<U>;
                fn $operator_fn(self, rhs: &Ss<U>) -> Self::Output {
                    let scalar_ss = Ss::<U>::new_from_scalar(self);
                    (&scalar_ss).$operator_fn(rhs)
                }
            }
        )*
    };
}

impl_scalar_math_operator_ss!([(Add, add), (Sub, sub), (Mul, mul), (Div, div)]);

macro_rules! impl_comb_ref_and_no_ref_operators {
    ([$(($operator:ident, $operator_fn:ident)), *]) => {
        $(
            impl<U: Time + 'static> $operator<&Ss<U>> for Ss<U>
            {
                type Output = Ss<U>;
                fn $operator_fn(self, rhs: &Ss<U>) -> Self::Output {
                    (&self).$operator_fn(rhs)
                }
            }

            impl<U: Time + 'static> $operator<Ss<U>> for &Ss<U>  
            {
                type Output = Ss<U>;
                fn $operator_fn(self, rhs: Ss<U>) -> Self::Output { 
                    self.$operator_fn(&rhs)
                }
            }
        )*
    };
}

impl_comb_ref_and_no_ref_operators!(
    [(Add, add), (Sub, sub), (Mul, mul), (Div, div)]
);

impl<U: Time + 'static> fmt::Display for Ss<U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\nA = {} B = {} C = {} D = {}",
            self.a(),
            self.b(),
            self.c(),
            self.d()
        )?;
        let time_domain = if TypeId::of::<U>() == TypeId::of::<Continuous>() {
            "Continuous"
        } else if TypeId::of::<U>() == TypeId::of::<Discrete>() {
            "Discrete"
        } else {
            "Unknown"
        };
        write!(f, "{time_domain}-time state-space model\n\n")
    }
}

#[cfg(test)]
mod tests {
    use core::f64;

    use approx::assert_abs_diff_eq;
    use num_complex::c64;
    use rand::Rng;

    use super::*;
    use crate::{
        FrequencyResponse,
        analysis::frequency_response::lin_space,
        systems::Tf,
        transformations::SsRealization::{ControllableCF, ObservableCF},
        utils::traits::Continuous,
    };

    fn rand_proper_tf<U: Rng>(
        rng: &mut U,
        max_order: usize,
    ) -> Tf<f64, Continuous> {
        let den_order = rng.random_range(1..=max_order);
        let num_order = rng.random_range(0..=den_order);

        let num: Vec<f64> = (0..=num_order)
            .map(|_| rng.random_range(-10.0..10.0))
            .collect();
        let mut den: Vec<f64> = (0..den_order)
            .map(|_| rng.random_range(-10.0..10.0))
            .collect();
        let mut den_max = rng.random_range(0.5..10.0); // ensure not too close to zero (possible division by close to zero)
        if rng.random_range(0..=1) != 1 {
            den_max *= -1.0; // negative for odd i
        }
        den.push(den_max);

        Tf::<f64, Continuous>::new(&num, &den)
    }

    #[test]
    fn new_ss() {
        let a = DMatrix::from_row_slice(2, 2, &[1., 0., 0., 1.]);
        let b = DMatrix::from_row_slice(2, 1, &[1., 1.]);
        let c = DMatrix::from_row_slice(1, 2, &[1., 0.]);
        let d = DMatrix::zeros(1, 1);

        assert!(
            Ss::<Continuous>::new(a.clone(), b.clone(), c.clone(), d.clone())
                .is_ok()
        );
        let a_wrong = DMatrix::zeros(1, 2);
        assert!(
            Ss::<Continuous>::new(
                a_wrong.clone(),
                b.clone(),
                c.clone(),
                d.clone()
            )
            .is_err()
        );
        let b_wrong = DMatrix::zeros(2, 2);
        assert!(
            Ss::<Continuous>::new(
                a.clone(),
                b_wrong.clone(),
                c.clone(),
                d.clone()
            )
            .is_err()
        );
        let c_wrong = DMatrix::zeros(2, 1);
        assert!(
            Ss::<Continuous>::new(
                a.clone(),
                b.clone(),
                c_wrong.clone(),
                d.clone()
            )
            .is_err()
        );
        let d_wrong = DMatrix::zeros(1, 2);
        assert!(
            Ss::<Continuous>::new(
                a.clone(),
                b.clone(),
                c.clone(),
                d_wrong.clone()
            )
            .is_err()
        );
    }

    #[test]
    fn ss_arithmetic() {
        let mut any_div_test = false;
        for _ in 0..1000 {
            let mut rng = rand::rng();
            let tf1 = rand_proper_tf(&mut rng, 5);
            let tf2 = rand_proper_tf(&mut rng, 5);

            let ss1 = tf1.to_ss_method(ControllableCF).unwrap();
            let ss2 = tf2.to_ss_method(ObservableCF).unwrap();

            let tf_add = tf1.clone() + tf2.clone();
            let tf_add = tf_add.to_ss().unwrap().to_tf().unwrap();
            let ss_add = ss1.clone() + ss2.clone();

            assert_abs_diff_eq!(
                tf_add,
                ss_add.to_tf().unwrap(),
                epsilon = 1e-2
            );

            let tf_sub = tf1.clone() - tf2.clone();
            let tf_sub = tf_sub.to_ss().unwrap().to_tf().unwrap();
            let ss_sub = ss1.clone() - ss2.clone();
            let ss_sub_tf = ss_sub.to_tf().unwrap();

            assert_abs_diff_eq!(tf_sub, ss_sub_tf, epsilon = 1e-2);

            let tf_mul = tf1.clone() * tf2.clone();
            let tf_mul = tf_mul.to_ss().unwrap().to_tf().unwrap();
            let ss_mul = ss1.clone() * ss2.clone();
            assert_abs_diff_eq!(
                tf_mul,
                ss_mul.to_tf().unwrap(),
                epsilon = 1e-2
            );

            if !tf2.is_strictly_proper() {
                any_div_test = true;
                let tf_div = tf1.clone() / tf2.clone();
                let tf_div = tf_div.to_ss().unwrap().to_tf().unwrap();
                let ss_div = ss1.clone() / ss2.clone();
                assert_abs_diff_eq!(
                    tf_div,
                    ss_div.to_tf().unwrap(),
                    epsilon = 1e-2
                );
            }
        }
        assert!(any_div_test);
    }

    #[test]
    fn ss_compund_assign() {
        let ss1 = (1.0 / Tf::s()).to_ss().unwrap();
        let ss2 = (Tf::s() / (1.0 + Tf::s())).to_ss().unwrap();

        let ss1_sum = ss1.clone() + ss2.clone();
        let mut ss1_sum_c = ss1.clone();
        ss1_sum_c += ss2.clone();
        assert_eq!(ss1_sum, ss1_sum_c);

        let ss1_sub = ss1.clone() - ss2.clone();
        let mut ss1_sub_c = ss1.clone();
        ss1_sub_c -= ss2.clone();
        assert_eq!(ss1_sub, ss1_sub_c);

        let ss1_mul = ss1.clone() * ss2.clone();
        let mut ss1_mul_c = ss1.clone();
        ss1_mul_c *= ss2.clone();
        assert_eq!(ss1_mul, ss1_mul_c);

        let ss1_div = ss1.clone() / ss2.clone();
        let mut ss1_div_c = ss1.clone();
        ss1_div_c /= ss2.clone();
        assert_eq!(ss1_div, ss1_div_c);
    }

    #[test]
    fn ss_scalar_arithmetic() {
        for scalar in lin_space(-10., 10., 1000) {
            if scalar.abs() < 0.1 {
                continue; // avoid division by zero
            }
            let tf_start = Tf::s() / (Tf::s() + 1.);
            let ss_start = tf_start.to_ss().unwrap();

            let ss_add = scalar + ss_start.clone();
            assert_abs_diff_eq!(
                ss_add.to_tf().unwrap(),
                scalar + tf_start.clone()
            );
            let ss_add = ss_start.clone() + scalar;
            assert_abs_diff_eq!(
                ss_add.to_tf().unwrap(),
                scalar + tf_start.clone()
            );

            let ss_sub = scalar - ss_start.clone();
            assert_abs_diff_eq!(
                ss_sub.to_tf().unwrap(),
                scalar - tf_start.clone()
            );
            let ss_sub = ss_start.clone() - scalar;
            assert_abs_diff_eq!(
                ss_sub.to_tf().unwrap(),
                tf_start.clone() - scalar
            );

            let ss_mul = scalar * ss_start.clone();
            assert_abs_diff_eq!(
                ss_mul.to_tf().unwrap(),
                scalar * tf_start.clone()
            );
            let ss_mul = ss_start.clone() * scalar;
            assert_abs_diff_eq!(
                ss_mul.to_tf().unwrap(),
                scalar * tf_start.clone()
            );

            let ss_div = scalar / ss_start.clone();
            assert_abs_diff_eq!(
                ss_div.to_tf().unwrap(),
                scalar / tf_start.clone()
            );
            let ss_div = ss_start.clone() / scalar;
            assert_abs_diff_eq!(
                ss_div.to_tf().unwrap(),
                tf_start.clone() / scalar
            );
        }
    }

    #[test]
    fn ss_connections() {
        let ss1 = (Tf::s() / (Tf::s() + 1.0)).to_ss().unwrap();
        let ss2 = (1.0 / Tf::s()).to_ss().unwrap();

        let ss_fb = ss1.clone().feedback(ss2.clone());
        let tf_fb = ss_fb.to_tf().unwrap();
        assert_abs_diff_eq!(tf_fb, Tf::s() / (Tf::s() + 2.0));

        let ss2 = Ss::<Continuous>::new_from_scalar(1.0);
        let ss_fb = ss1.clone().feedback(ss2.clone());
        let tf_fb = ss_fb.to_tf().unwrap();
        assert_abs_diff_eq!(tf_fb, 0.5 * Tf::s() / (Tf::s() + 0.5));

        let ss_add = ss1.clone() + ss2.clone();
        let ss_parallel = ss1.parallel(&ss2);
        assert_eq!(ss_add, ss_parallel);

        let ss_mul = ss2.clone() * ss1.clone();
        let ss_series = ss1.series(&ss2);
        assert_eq!(ss_mul, ss_series);
    }

    #[test]
    fn ss_and_tf_feedback_match() {
        for _ in 0..1000 {
            let mut rng = rand::rng();

            let tf1 = rand_proper_tf(&mut rng, 5);
            let tf2 = rand_proper_tf(&mut rng, 5);

            let ss1 = tf1.to_ss_method(ControllableCF).unwrap();
            let ss2 = tf2.to_ss_method(ObservableCF).unwrap();

            let ss_fb = ss1.feedback(ss2);
            let ss_fb = ss_fb.to_tf().unwrap();
            let tf_fb = tf1.feedback(&tf2);

            assert_abs_diff_eq!(tf_fb, ss_fb, epsilon = 1e-1);
        }
    }

    #[test]
    fn ss_display() {
        let sys = (1.0 / Tf::s()).to_ss().unwrap();
        println!("1: {sys}");

        let sys = Ss::<Discrete>::new_from_scalar(2.0);
        println!("2: {sys}");

        let sys = ((1.0 + Tf::s()) / ((Tf::s() - 1.0) * Tf::s()).powi(2))
            .to_ss()
            .unwrap();
        println!("3: {sys}after");
    }

    #[test]
    fn ss_reference_arithmetic() {
        let mut ss = (1.0 / Tf::s().powi(2)).to_ss().unwrap();
        let ss_org = ss.clone();

        ss += &ss_org;
        ss += &ss_org;

        let ans = (ss_org.clone() + ss_org.clone()) + ss_org.clone();

        let freq = c64(0.0, 1.2);
        let resp1 = ss.freq_response(&freq);
        let resp2 = ans.freq_response(&freq);

        assert_abs_diff_eq!(resp1.re, resp2.re, epsilon = 1e-9);
        assert_abs_diff_eq!(resp1.im, resp2.im, epsilon = 1e-9);


        let s_inv = (1.0 / Tf::s()).to_ss().unwrap();
        let _ = 1.0 + &s_inv * &s_inv / (1.0 + &s_inv / 1.0);
    }
}
