use std::{error::Error, marker::PhantomData};

use crate::{
    slicotrs::ss2tf_tb04ad,
    tf::Tf,
    traits::{Time, Zero},
};
extern crate nalgebra as na;
use na::DMatrix;

#[derive(Debug, Clone)]
pub struct Ss<U: Time> {
    a: DMatrix<f64>,
    b: DMatrix<f64>,
    c: DMatrix<f64>,
    d: DMatrix<f64>,
    time: PhantomData<U>,
}

impl<U: Time + 'static> Ss<U> {
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

    pub fn a(&self) -> &DMatrix<f64> {
        &self.a
    }

    pub fn b(&self) -> &DMatrix<f64> {
        &self.b
    }

    pub fn c(&self) -> &DMatrix<f64> {
        &self.c
    }

    pub fn d(&self) -> &DMatrix<f64> {
        &self.d
    }

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
