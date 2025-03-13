use crate::{
    slicotrs::ss2tf_tb04ad,
    ss::Ss,
    tf::Tf,
    traits::{Time, Zero},
};
use nalgebra::DMatrix;
use std::error::Error;

pub fn ss2tf<U: Time + 'static>(
    ss: &Ss<U>,
) -> Result<Tf<f64, U>, Box<dyn Error + 'static>> {
    ss2tf_mat(ss.a(), ss.b(), ss.c(), ss.d())
}

pub fn ss2tf_mat<U: Time + 'static>(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
) -> Result<Tf<f64, U>, Box<dyn Error + 'static>> {
    ss2tf_tb04ad::<U>(a, b, c, d)
}

pub fn tf2ss<U: Time + 'static>(
    tf: Tf<f64, U>,
    method: SsRealization,
) -> Result<Ss<U>, Box<dyn Error + 'static>> {
    match method {
        SsRealization::ControllableCF => tf2ss_controllable(tf),
        SsRealization::ObservableCF => tf2ss_observable(tf),
    }
}

fn tf2ss_observable<U: Time + 'static>(
    tf: Tf<f64, U>,
) -> Result<Ss<U>, Box<dyn Error + 'static>> {
    if !tf.is_proper() {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Transfer function must be proper",
        )));
    }

    let tf = tf.normalize();

    let (_, den_deg) = tf.degree_num_den();
    let nx = den_deg;

    let mut a = DMatrix::zeros(nx, nx);
    a.view_mut((1, 0), (nx - 1, nx - 1))
        .copy_from(&DMatrix::identity(nx - 1, nx - 1));
    for row in 0..nx {
        a[(row, nx - 1)] = -tf.denominator()[row];
    }

    let mut d = DMatrix::zeros(1, 1);
    if !tf.is_strictly_proper() {
        assert!(tf.numerator().len() > nx);
        d[(0, 0)] = tf.numerator()[nx];
    }

    let mut num_extended = tf.numerator().to_vec();
    num_extended.resize(nx, 0.);

    let b_values: Vec<f64> = if d[(0, 0)].is_zero() {
        num_extended
    } else {
        let d_const = d[(0, 0)];
        let num_vals = &num_extended;
        let den_vals = &tf.numerator()[..nx];
        num_vals
            .iter()
            .zip(den_vals.iter())
            .map(|(num_i, den_i)| num_i - d_const * den_i)
            .collect()
    };

    let b = DMatrix::from_column_slice(nx, 1, &b_values);

    let mut c = DMatrix::zeros(1, nx);
    c[(0, 0)] = 1.;

    Ok(Ss::<U>::new(a, b, c, d)?)
}

fn tf2ss_controllable<U: Time + 'static>(
    tf: Tf<f64, U>,
) -> Result<Ss<U>, Box<dyn Error + 'static>> {
    let ss = tf2ss_observable(tf)?;

    let a = ss.a().transpose();
    let b = ss.c().transpose();
    let c = ss.b().transpose();
    let d = ss.d().transpose();

    Ok(Ss::<U>::new(a, b, c, d)?)
}

#[derive(Debug)]
pub enum SsRealization {
    ObservableCF,
    ControllableCF,
}

impl Default for SsRealization {
    fn default() -> Self {
        SsRealization::ObservableCF
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::traits::Continuous;
    use approx::assert_abs_diff_eq;
    #[test]
    fn ss2tf_test() {
        let a = DMatrix::from_row_slice(2, 2, &[0., 1., 0., 0.]);
        let b = DMatrix::from_row_slice(2, 1, &[0., 1.]);
        let c = DMatrix::from_row_slice(1, 2, &[1., 0.]);
        let d = DMatrix::zeros(1, 1);

        let ss = Ss::<Continuous>::new(a, b, c, d).unwrap();
        let tf = ss2tf(&ss).unwrap();

        let tf_ans = 1./Tf::s().powi(2);
        println!("ss2Tf: \n{}", tf);
        assert_abs_diff_eq!(tf, tf_ans);
    }

    #[test]

    fn tf2ss_test() {
        let tf = 1./Tf::s().powi(2);

        let ss = tf2ss(tf, SsRealization::ObservableCF).unwrap();
        println!("{:?}", ss);
        let tf_ret = ss2tf(&ss).unwrap();
        println!("tf ret: {}\n", tf_ret);
        let a_ans = DMatrix::from_row_slice(2, 2, &[0., 1., 0., 0.]);
        let b_ans = DMatrix::from_row_slice(2, 1, &[0., 1.]);
        let c_ans = DMatrix::from_row_slice(1, 2, &[1., 0.]);
        let d_ans = DMatrix::from_row_slice(1, 1, &[0.]);
        assert_abs_diff_eq!(ss.a(), &a_ans);
        assert_abs_diff_eq!(ss.b(), &b_ans);
        assert_abs_diff_eq!(ss.c(), &c_ans);
        assert_abs_diff_eq!(ss.d(), &d_ans);

    }
}
