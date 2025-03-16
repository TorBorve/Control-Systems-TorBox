use std::{error::Error, marker::PhantomData};

use crate::traits::Time;
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
#[derive(Debug, Clone)]
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
        Ss::<U>::verify_dimensions(&a, &b, &c, &d)?;
        Ok(Self {
            a,
            b,
            c,
            d,
            time: PhantomData::<U>,
        })
    }

    /// Returns a reference to the state matrix (A).
    pub fn a(&self) -> &DMatrix<f64> {
        &self.a
    }

    /// Returns a reference to the state matrix (B).
    pub fn b(&self) -> &DMatrix<f64> {
        &self.b
    }

    /// Returns a reference to the state matrix (C).
    pub fn c(&self) -> &DMatrix<f64> {
        &self.c
    }

    /// Returns a reference to the state matrix (D).
    pub fn d(&self) -> &DMatrix<f64> {
        &self.d
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
}

#[cfg(test)]
mod tests {
    use crate::traits::Continuous;

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
}
