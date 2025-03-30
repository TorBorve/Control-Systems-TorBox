use std::{
    error::Error,
    marker::PhantomData,
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use crate::{
    slicotrs::{h2_norm, hinf_norm, zeros},
    traits::Time,
};
extern crate nalgebra as na;
use na::DMatrix;
use num_complex::Complex64;

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
            "Ss is invalid: Message: {}, Ss: {:?}",
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
    pub fn series(self, sys2: Self) -> Self {
        assert_eq!(
            self.noutputs(),
            sys2.ninputs(),
            "Self output must be same as sys2 inputs"
        );
        self.assert_valid();
        sys2.assert_valid();
        let ret = sys2 * self;
        ret.assert_valid();
        ret
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
        assert_eq!(
            self.shape(),
            sys2.shape(),
            "System must have same shape to connnect in parallel"
        );
        self.assert_valid();
        sys2.assert_valid();
        let ret = self + sys2;
        ret.assert_valid();
        ret
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

    /// Computes the H2 norm of the system.
    ///
    /// The H2 norm is a measure of the energy gain from the input to the output
    /// in the frequency domain. It is defined as:
    ///
    /// ```txt
    /// ||G||_H2 = sqrt( trace( C * W_c * C^T ) )
    /// ```
    ///
    /// where `W_c` is the controllability Gramian.
    ///
    /// # Returns
    /// - `Ok(f64)`: The H2 norm of the system.
    /// - `Err(String)`: An error message if the computation fails.
    pub fn norm_h2(self) -> Result<f64, String> {
        h2_norm::<U>(self.a, self.b, self.c, self.d).map_err(|e| e.to_string())
    }

    /// Computes the H∞ (H-infinity) norm of the system.
    ///
    /// The H∞ norm represents the maximum singular value of the transfer
    /// function across all frequencies. It quantifies the worst-case
    /// amplification of the system from input to output.
    ///
    /// # Returns
    /// - `Ok(f64)`: The H∞ norm of the system.
    /// - `Err(String)`: An error message if the computation fails.
    pub fn norm_hinf(self) -> Result<f64, String> {
        hinf_norm::<U>(self.a, self.b, self.c, self.d)
    }

    /// Computes the poles of the system
    ///
    /// The poles of the system are equal to the eigen values of A.
    ///
    /// # Returns
    /// - Vector with complex eigen values
    pub fn poles(&self) -> Vec<Complex64> {
        let eigen_values = self.a().complex_eigenvalues();
        assert_eq!(eigen_values.ncols(), 1);
        eigen_values.as_slice().to_vec()
    }

    /// Computes the invariant zeros of a continuous-time state-space system.
    ///
    /// # Returns
    /// A `Result` containing a vector of complex numbers representing the
    /// system zeros, or an error message.
    pub fn zeros(&self) -> Result<Vec<Complex64>, String> {
        zeros(self.a(), self.b(), self.c(), self.d())
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
        assert_eq!(self.shape(), rhs.shape());
        self.assert_valid();
        rhs.assert_valid();

        let nu = self.ninputs();
        let ny = self.noutputs();

        let nx1 = self.order();
        let nx2 = rhs.order();

        let mut a_new = DMatrix::zeros(nx1 + nx2, nx1 + nx2);
        a_new.view_mut((0, 0), (nx1, nx1)).copy_from(self.a());
        a_new.view_mut((nx1, nx1), (nx2, nx2)).copy_from(rhs.a());

        let mut b_new = DMatrix::zeros(nx1 + nx2, nu);
        b_new.view_mut((0, 0), (nx1, nu)).copy_from(self.b());
        b_new.view_mut((nx1, 0), (nx2, nu)).copy_from(rhs.b());

        let mut c_new = DMatrix::zeros(ny, nx1 + nx2);
        c_new.view_mut((0, 0), (ny, nx1)).copy_from(self.c());

        c_new.view_mut((0, nx1), (ny, nx2)).copy_from(rhs.c());
        let d_new = self.d() + rhs.d();

        Ss::<U>::new(a_new, b_new, c_new, d_new).unwrap()
    }
}

impl<U: Time + 'static> Neg for Ss<U> {
    type Output = Ss<U>;

    /// Negates the output of a state-space system.
    ///
    /// I.e. y_new = -y    
    fn neg(mut self) -> Self::Output {
        self.c = -self.c();
        self.d = -self.d();
        self.assert_valid();
        self
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
        self.assert_valid();
        rhs.assert_valid();

        let nx2 = self.order();
        let nu2 = self.ninputs();
        let ny2 = self.noutputs();
        let nx1 = rhs.order();
        let nu1 = rhs.ninputs();
        let ny1 = rhs.noutputs();
        assert_eq!(nu2, ny1);

        let mut a_new = DMatrix::zeros(nx1 + nx2, nx1 + nx2);
        a_new.view_mut((0, 0), (nx1, nx1)).copy_from(rhs.a());
        a_new
            .view_mut((nx1, 0), (nx2, nx1))
            .copy_from(&(self.b() * rhs.c()));
        a_new.view_mut((nx1, nx1), (nx2, nx2)).copy_from(self.a());

        let mut b_new = DMatrix::zeros(nx1 + nx2, nu1);
        b_new.view_mut((0, 0), (nx1, nu1)).copy_from(rhs.b());
        b_new
            .view_mut((nx1, 0), (nx2, nu1))
            .copy_from(&(self.b() * rhs.d()));

        let mut c_new = DMatrix::zeros(ny2, nx1 + nx2);
        c_new
            .view_mut((0, 0), (ny2, nx1))
            .copy_from(&(self.d() * rhs.c()));
        c_new.view_mut((0, nx1), (ny2, nx2)).copy_from(self.c());

        let d_new = self.d() * rhs.d();

        Ss::<U>::new(a_new, b_new, c_new, d_new).unwrap()
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
        self.assert_valid();
        rhs.assert_valid();
        self * rhs.inv().unwrap()
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
        )*
    };
}

impl_scalar_math_operator_ss!([(Add, add), (Sub, sub), (Mul, mul), (Div, div)]);

#[cfg(test)]
mod tests {
    use core::f64;

    use approx::assert_abs_diff_eq;
    use num_complex::c64;
    use rand::Rng;

    use crate::{
        Continuous, Tf, lin_space, ss2tf,
        tests::rand_proper_tf,
        tf2ss,
        transforms::SsRealization::{ControllableCF, ObservableCF},
    };

    use super::*;

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

            let ss1 = tf2ss(tf1.clone(), ControllableCF).unwrap();
            let ss2 = tf2ss(tf2.clone(), ObservableCF).unwrap();

            let tf_add = tf1.clone() + tf2.clone();
            let tf_add =
                ss2tf(&tf2ss(tf_add, ControllableCF).unwrap()).unwrap();
            let ss_add = ss1.clone() + ss2.clone();

            assert_abs_diff_eq!(
                tf_add,
                ss2tf(&ss_add).unwrap(),
                epsilon = 1e-2
            );

            let tf_sub = tf1.clone() - tf2.clone();
            let tf_sub =
                ss2tf(&tf2ss(tf_sub, ControllableCF).unwrap()).unwrap();
            let ss_sub = ss1.clone() - ss2.clone();
            let ss_sub_tf = ss2tf(&ss_sub).unwrap();

            assert_abs_diff_eq!(tf_sub, ss_sub_tf, epsilon = 1e-2);

            let tf_mul = tf1.clone() * tf2.clone();
            let tf_mul = ss2tf(&tf2ss(tf_mul, ObservableCF).unwrap()).unwrap();
            let ss_mul = ss1.clone() * ss2.clone();
            assert_abs_diff_eq!(
                tf_mul,
                ss2tf(&ss_mul).unwrap(),
                epsilon = 1e-2
            );

            if !tf2.is_strictly_proper() {
                any_div_test = true;
                let tf_div = tf1.clone() / tf2.clone();
                let tf_div =
                    ss2tf(&tf2ss(tf_div, ObservableCF).unwrap()).unwrap();
                let ss_div = ss1.clone() / ss2.clone();
                assert_abs_diff_eq!(
                    tf_div,
                    ss2tf(&ss_div).unwrap(),
                    epsilon = 1e-2
                );
            }
        }
        assert!(any_div_test);
    }

    #[test]
    fn ss_compund_assign() {
        let ss1 = tf2ss(1.0 / Tf::s(), ObservableCF).unwrap();
        let ss2 = tf2ss(Tf::s() / (1.0 + Tf::s()), ObservableCF).unwrap();

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
            let ss_start = tf2ss(tf_start.clone(), ObservableCF).unwrap();

            let ss_add = scalar + ss_start.clone();
            assert_abs_diff_eq!(
                ss2tf(&ss_add).unwrap(),
                scalar + tf_start.clone()
            );
            let ss_add = ss_start.clone() + scalar;
            assert_abs_diff_eq!(
                ss2tf(&ss_add).unwrap(),
                scalar + tf_start.clone()
            );

            let ss_sub = scalar - ss_start.clone();
            assert_abs_diff_eq!(
                ss2tf(&ss_sub).unwrap(),
                scalar - tf_start.clone()
            );
            let ss_sub = ss_start.clone() - scalar;
            assert_abs_diff_eq!(
                ss2tf(&ss_sub).unwrap(),
                tf_start.clone() - scalar
            );

            let ss_mul = scalar * ss_start.clone();
            assert_abs_diff_eq!(
                ss2tf(&ss_mul).unwrap(),
                scalar * tf_start.clone()
            );
            let ss_mul = ss_start.clone() * scalar;
            assert_abs_diff_eq!(
                ss2tf(&ss_mul).unwrap(),
                scalar * tf_start.clone()
            );

            let ss_div = scalar / ss_start.clone();
            assert_abs_diff_eq!(
                ss2tf(&ss_div).unwrap(),
                scalar / tf_start.clone()
            );
            let ss_div = ss_start.clone() / scalar;
            assert_abs_diff_eq!(
                ss2tf(&ss_div).unwrap(),
                tf_start.clone() / scalar
            );
        }
    }

    #[test]
    fn ss_connections() {
        let ss1 = tf2ss(Tf::s() / (Tf::s() + 1.0), ObservableCF).unwrap();
        let ss2 = tf2ss(1.0 / Tf::s(), ObservableCF).unwrap();

        let ss_fb = ss1.clone().feedback(ss2.clone());
        let tf_fb = ss2tf(&ss_fb).unwrap();
        assert_abs_diff_eq!(tf_fb, Tf::s() / (Tf::s() + 2.0));

        let ss2 = Ss::<Continuous>::new_from_scalar(1.0);
        let ss_fb = ss1.clone().feedback(ss2.clone());
        let tf_fb = ss2tf(&ss_fb).unwrap();
        assert_abs_diff_eq!(tf_fb, 0.5 * Tf::s() / (Tf::s() + 0.5));

        let ss_add = ss1.clone() + ss2.clone();
        let ss_parallel = ss1.clone().parallel(ss2.clone());
        assert_eq!(ss_add, ss_parallel);

        let ss_mul = ss2.clone() * ss1.clone();
        let ss_series = ss1.series(ss2);
        assert_eq!(ss_mul, ss_series);
    }

    #[test]
    fn ss_and_tf_feedback_match() {
        for _ in 0..1000 {
            let mut rng = rand::rng();

            let tf1 = rand_proper_tf(&mut rng, 5);
            let tf2 = rand_proper_tf(&mut rng, 5);

            let ss1 = tf2ss(tf1.clone(), ObservableCF).unwrap();
            let ss2 = tf2ss(tf2.clone(), ObservableCF).unwrap();

            let ss_fb = ss1.feedback(ss2);
            let ss_fb = ss2tf(&ss_fb).unwrap();
            let tf_fb = tf1.feedback(tf2);

            assert_abs_diff_eq!(tf_fb, ss_fb, epsilon = 1e-1);
        }
    }

    #[test]
    fn ss_system_norms() {
        let sys_tf = Tf::s() / (Tf::s() + 1.0);
        let sys = tf2ss(sys_tf, ObservableCF).unwrap();
        assert!(sys.clone().norm_h2().is_err());
        assert_abs_diff_eq!(sys.clone().norm_hinf().unwrap(), 1.0);

        let sys = 2.0 / (Tf::s() + 1.0);
        let sys = tf2ss(sys, ObservableCF).unwrap();
        assert_abs_diff_eq!(
            sys.clone().norm_h2().unwrap(),
            2.0 / 2.0_f64.sqrt()
        );
        assert_abs_diff_eq!(sys.clone().norm_hinf().unwrap(), 2.0);

        let damps = [0.1, 0.2, 0.3, 0.4];
        let matlab_results = [5.0252, 2.5515, 1.7471, 1.3639];
        for (expected_hinf_norm, damping) in
            matlab_results.iter().zip(damps.iter())
        {
            let sys_tf =
                1.0 / (Tf::s().powi(2) + 2.0 * damping * Tf::s() + 1.0);
            let sys = tf2ss(sys_tf.clone(), ObservableCF).unwrap();
            assert_abs_diff_eq!(
                sys.norm_hinf().unwrap(),
                expected_hinf_norm,
                epsilon = 1e-2
            );
        }
    }

    #[test]
    fn ss_poles() {
        let mut rng = rand::rng();
        for _ in 0..10 {
            let mut poles = vec![];
            let mut sys_tf = Tf::new_from_scalar(1.0);
            for new_pole in (0..3).map(|_| {
                c64(rng.random_range(-10.0..10.0), rng.random_range(0.0..10.0))
            }) {
                if new_pole.norm() < 1e-1 {
                    continue;
                }
                let new_pole_conj = new_pole.conj();
                poles.push(new_pole);
                poles.push(new_pole_conj);

                sys_tf *= 1.0
                    / (Tf::s().powi(2) - 2.0 * new_pole.re * Tf::s()
                        + (new_pole.im.powi(2) + new_pole.re.powi(2)));
                let calc_poles =
                    tf2ss(sys_tf.clone(), ObservableCF).unwrap().poles();

                assert_eq!(calc_poles.len(), poles.len());
                for pole in &poles {
                    let min_distance = calc_poles
                        .iter()
                        .map(|p| (p - pole).norm())
                        .fold(f64::INFINITY, f64::min);
                    assert_abs_diff_eq!(min_distance, 0.0, epsilon = 1e-1);
                }
            }
        }

        let ss = tf2ss(1.0 / Tf::s(), ObservableCF).unwrap();
        let poles = ss.poles();
        assert_eq!(poles.len(), 1);
        let pole = poles[0];
        assert_abs_diff_eq!(pole.re, 0.0);
        assert_abs_diff_eq!(pole.im, 0.0);
    }

    #[test]
    fn ss_zeros() {
        let tf = (Tf::s() - 1.0) * (Tf::s() + 4.0) / (Tf::s() + 2.0).powi(2);
        println!("tf: \n{}", tf);
        let zeros = tf2ss(tf, ControllableCF).unwrap().zeros().unwrap();
        println!("zeros: {:?}", zeros);
        assert_eq!(zeros.len(), 2);

        let tf = 1.0 / Tf::s();
        let zeros = tf2ss(tf, ObservableCF).unwrap().zeros().unwrap();
        assert_eq!(zeros.len(), 0);

        for _ in 0..100 {
            let mut rng = rand::rng();
            for zero in (0..100).map(|_| rng.random_range(-100.0..100.0)) {
                let tf = (Tf::s() - zero) / (Tf::s() + 1.0).powi(2);
                let sys = tf2ss(tf, ObservableCF).unwrap();
                let zeros = sys.zeros().unwrap();
                assert_eq!(zeros.len(), 1);
                assert_abs_diff_eq!(zeros[0].re, zero, epsilon = 1e-2);
                assert_abs_diff_eq!(zeros[0].im, 0.0, epsilon = 1e-2);
            }

            let mut tf = Tf::new_from_scalar(1.0);
            let mut zeros = vec![];
            for new_zero in (0..3).map(|_| {
                c64(rng.random_range(-10.0..-1.0), rng.random_range(0.0..10.0))
            }) {
                if new_zero.norm() < 1e-1 {
                    continue;
                }
                let new_zero = c64(new_zero.re, 0.0);
                zeros.push(new_zero);
                zeros.push(new_zero.conj());
                tf *= (Tf::s().powi(2) - 2.0 * new_zero.re * Tf::s()
                    + new_zero.norm_sqr())
                    / (Tf::s().powi(2) + 0.0);
                let sys = tf2ss(tf.clone(), ObservableCF).unwrap();
                let calc_zeros = sys.zeros().unwrap();
                assert_eq!(zeros.len(), calc_zeros.len());
                for zero in &zeros {
                    let min_dist = calc_zeros
                        .iter()
                        .map(|z| (z - zero).norm())
                        .fold(f64::INFINITY, f64::min);
                    assert_abs_diff_eq!(min_dist, 0.0, epsilon = 1e-1);
                }
            }
        }
    }
}
