use crate::{slicotrs::ss2tf_tb04ad, ss::Ss, tf::Tf, traits::Time};
use nalgebra::DMatrix;
use std::error::Error;

/// Converts a state-space representation to a transfer function.
///
/// # Arguments
///
/// * `ss` - A reference to a state-space system.
///
/// # Returns
///
/// Returns a `Result` containing a transfer function (`Tf<f64, U>`) or an error
/// if the conversion fails.
///
/// # Example
///
/// ```rust
/// use control_systems_torbox::*;
/// let ss = tf2ss(1.0/Tf::s(), SsRealization::ControllableCF).unwrap();
/// let tf = ss2tf::<Continuous>(&ss);
/// ```
pub fn ss2tf<U: Time + 'static>(
    ss: &Ss<U>,
) -> Result<Tf<f64, U>, Box<dyn Error + 'static>> {
    ss2tf_mat(ss.a(), ss.b(), ss.c(), ss.d())
}

/// Converts state-space matrices to a transfer function.
///
/// # Arguments
///
/// * `a` - State matrix.
/// * `b` - Input matrix.
/// * `c` - Output matrix.
/// * `d` - Feedthrough matrix.
///
/// # Returns
///
/// Returns a `Result` containing a transfer function (`Tf<f64, U>`) or an error
/// if the conversion fails.
///
/// # Example
///
/// ```rust
/// use control_systems_torbox::*;
/// use nalgebra::DMatrix;
/// let a = DMatrix::identity(1, 1);
/// let b = DMatrix::identity(1, 1);
/// let c = DMatrix::identity(1, 1);
/// let d = DMatrix::zeros(1, 1);
/// let tf = ss2tf_mat::<Discrete>(&a, &b, &c, &d);
/// ```
pub fn ss2tf_mat<U: Time + 'static>(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
) -> Result<Tf<f64, U>, Box<dyn Error + 'static>> {
    ss2tf_tb04ad::<U>(a, b, c, d)
}

/// Converts a transfer function to a state-space representation.
///
/// # Arguments
///
/// * `tf` - A transfer function (`Tf<f64, U>`).
/// * `method` - The realization method used for the conversion such as
///   Controllable Canonical Form or Observable Canonical Form.
///
/// # Returns
///
/// Returns a `Result` containing a state-space system (`Ss<U>`) or an error if
/// the conversion fails.
///
/// # Example
///
/// ```rust
/// use control_systems_torbox::*;
/// let tf = 1.0/Tf::s();
/// let ss = tf2ss(tf, SsRealization::ControllableCF);
/// ```
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

    let mut b_values = num_extended;
    assert!(tf.denominator().len() - 1 <= b_values.len());
    for (i, den_i) in tf.denominator().iter().enumerate().take(nx) {
        b_values[i] += -d[(0, 0)] * den_i;
    }

    let b = DMatrix::from_column_slice(nx, 1, &b_values);

    let mut c = DMatrix::zeros(1, nx);
    c[(0, nx - 1)] = 1.;

    Ss::<U>::new(a, b, c, d)
}

fn tf2ss_controllable<U: Time + 'static>(
    tf: Tf<f64, U>,
) -> Result<Ss<U>, Box<dyn Error + 'static>> {
    let ss = tf2ss_observable(tf)?;

    let a = ss.a().transpose();
    let b = ss.c().transpose();
    let c = ss.b().transpose();
    let d = ss.d().transpose();

    Ss::<U>::new(a, b, c, d)
}

/// Enumeration for the different realizations of state-space systems.
///
/// # Variants
///
/// * `ObservableCF` - Observable canonical form.
/// * `ControllableCF` - Controllable canonical form.
#[derive(Debug, Clone, Copy, Default)]
pub enum SsRealization {
    #[default]
    ObservableCF,
    ControllableCF,
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use rayon::prelude::*;

    use crate::traits::Continuous;
    use approx::assert_abs_diff_eq;
    use rand::{Rng, seq::IteratorRandom};
    #[test]
    fn ss2tf_test() {
        let a = DMatrix::from_row_slice(2, 2, &[0., 1., 0., 0.]);
        let b = DMatrix::from_row_slice(2, 1, &[0., 1.]);
        let c = DMatrix::from_row_slice(1, 2, &[1., 0.]);
        let d = DMatrix::zeros(1, 1);

        let ss = Ss::<Continuous>::new(a, b, c, d).unwrap();
        let tf = ss2tf(&ss).unwrap();

        let tf_ans = 1. / Tf::s().powi(2);
        println!("ss2Tf: \n{}", tf);
        assert_abs_diff_eq!(tf, tf_ans);
    }

    #[test]

    fn tf2ss_test() {
        let tf = 1. / Tf::s().powi(2);

        let ss = tf2ss(tf.clone(), SsRealization::ObservableCF).unwrap();
        let tf_ret = ss2tf(&ss).unwrap();
        let a_ans = DMatrix::from_row_slice(2, 2, &[0., 0., 1., 0.]);
        let b_ans = DMatrix::from_row_slice(2, 1, &[1., 0.]);
        let c_ans = DMatrix::from_row_slice(1, 2, &[0., 1.]);
        let d_ans = DMatrix::from_row_slice(1, 1, &[0.]);
        assert_abs_diff_eq!(ss.a(), &a_ans);
        assert_abs_diff_eq!(ss.b(), &b_ans);
        assert_abs_diff_eq!(ss.c(), &c_ans);
        assert_abs_diff_eq!(ss.d(), &d_ans);

        assert_abs_diff_eq!(tf_ret.normalize(), tf.normalize());
    }

    #[test]
    fn convert_between_tf_and_ss() {
        let start = Instant::now();
        let num_iter = 10000;
        (0..num_iter).into_par_iter().for_each(|i| {
            let mut rng = rand::rng();
            let den_order = rng.random_range(1..=20) as usize;
            let num_order = rng.random_range(0..=den_order);

            let num: Vec<f64> = (0..=num_order)
                .map(|_| rng.random_range(-10.0..10.0))
                .collect();
            let mut den: Vec<f64> = (0..=den_order - 1)
                .map(|_| rng.random_range(-10.0..10.0))
                .collect();
            let mut den_max = rng.random_range(0.5..10.0); // ensure not too close to zero (possible division by close to zero)
            if i % 2 != 0 {
                den_max *= -1.0; // negative for odd i
            }
            den.push(den_max);

            let tf = Tf::<f64, Continuous>::new(&num, &den);
            let methods =
                [SsRealization::ControllableCF, SsRealization::ObservableCF];
            let ss =
                tf2ss(tf.clone(), *methods.iter().choose(&mut rng).unwrap())
                    .unwrap();
            let tf_ret = ss2tf(&ss).unwrap().normalize();

            let tf = tf.normalize();
            assert_abs_diff_eq!(tf, tf_ret, epsilon = 1e-3);
        });
        println!(
            "Time ss2tf(tf2ss()) transforms avg: {:?}",
            start.elapsed() / num_iter
        );
    }
}
