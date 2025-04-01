extern crate libc;
extern crate netlib_src; // enable linking with blas and lapack
// use crate::{Continuous, Discrete, ss::Ss, tf::Tf, traits::Time};
use crate::{
    slicot_wrapper::{ab08nd_, ab13bd_, ab13dd_, dggev_},
    systems::Ss,
    utils::traits::{Continuous, Discrete, Time},
};
use libc::{c_double, c_int};
use nalgebra::DMatrix;
use num_complex::{Complex64, c64};
use std::{error::Error, ffi::CString};

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
/// - `Err(dyn Error)`: An error message if the computation fails.
fn h2_norm<U: Time + 'static>(
    mut a: DMatrix<f64>,
    mut b: DMatrix<f64>,
    mut c: DMatrix<f64>,
    mut d: DMatrix<f64>,
) -> Result<f64, Box<dyn Error + 'static>> {
    let time_domain = if std::any::TypeId::of::<U>()
        == std::any::TypeId::of::<Continuous>()
    {
        CString::new("C").unwrap()
    } else if std::any::TypeId::of::<U>() == std::any::TypeId::of::<Discrete>()
    {
        CString::new("D").unwrap()
    } else {
        unreachable!("Time domain must be Continuous or Discrete");
    };

    let h2_or_l2 = CString::new("H").unwrap();

    Ss::<U>::verify_dimensions(&a, &b, &c, &d)?;

    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();

    let lda = a.nrows();
    let ldb = b.nrows();
    let ldc = c.nrows();
    let ldd = d.nrows();

    let mut nq = -1 as c_int;
    let tol = -1 as c_double;

    let ldwork = 5 * 1
        .max(m * (n + m) + (n * (n + 5)).max((m * (m + 2)).max(4 * p)))
        .max(n * (n.max(p) + 4) + n.min(p));
    let mut dwork = vec![0.0 as c_double; ldwork];

    let mut iwarn = -1 as c_int;
    let mut info = -1 as c_int;

    let h2_norm = unsafe {
        ab13bd_(
            time_domain.as_ptr(),
            h2_or_l2.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            a.as_mut_ptr(),
            &(lda as c_int),
            b.as_mut_ptr(),
            &(ldb as c_int),
            c.as_mut_ptr(),
            &(ldc as c_int),
            d.as_mut_ptr(),
            &(ldd as c_int),
            &mut nq,
            &tol,
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            &mut iwarn,
            &mut info,
        )
    };

    if info != 0 {
        return Err(Box::new(std::io::Error::other(format!(
            "SLICOT ab13bd failed with the follwing error indicator: {}",
            info
        ))));
    } else if iwarn != 0 {
        return Err(Box::new(std::io::Error::other(format!(
            "SLICOT ab13bd returned with warning about numerical stability. Iwarn = {}",
            iwarn
        ))));
    }

    Ok(h2_norm)
}

/// Computes the H∞ (H-infinity) norm of the system.
///
/// The H∞ norm represents the maximum singular value of the transfer function
/// across all frequencies. It quantifies the worst-case amplification of the
/// system from input to output.
///
/// # Returns
/// - `Ok(f64)`: The H∞ norm of the system.
/// - `Err(String)`: An error message if the computation fails.
fn hinf_norm<U: Time + 'static>(
    mut a: DMatrix<f64>,
    mut b: DMatrix<f64>,
    mut c: DMatrix<f64>,
    mut d: DMatrix<f64>,
) -> Result<f64, String> {
    Ss::<U>::verify_dimensions(&a, &b, &c, &d).map_err(|e| e.to_string())?;

    let time_domain = if std::any::TypeId::of::<U>()
        == std::any::TypeId::of::<Continuous>()
    {
        CString::new("C").unwrap()
    } else if std::any::TypeId::of::<U>() == std::any::TypeId::of::<Discrete>()
    {
        CString::new("D").unwrap()
    } else {
        unreachable!("Time domain must be Continuous or Discrete");
    };

    let e_shape = CString::new("I").unwrap(); // E is identity
    let equilibration = CString::new("S").unwrap(); // Perform scaling
    let d_is_nonzero = CString::new("D").unwrap();

    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();

    let mut e = DMatrix::identity(n, n);

    let mut f_peak = [0.0, 1.0]; // initial guess not active
    let mut peak_gain = [0.0, 0.0];
    let tol = 1e-6;

    let mut iwork = vec![0; n];
    let ldwork = 10
        * 1.max(
            15 * n * n
                + p * p
                + m * m
                + (6 * n + 3) * (p + m)
                + 4 * p * m
                + n * m
                + 22 * n
                + 7 * p.min(m),
        );
    let mut dwork = vec![0.0; ldwork];

    let lcwork = 10 * 1.max((n + m) * (n + p) + 2 * p.min(m) + p.max(m));
    let mut cwork = vec![0.; 2 * lcwork]; // Mulipy by two because cwork is complex 16 bytes not f64 // cwork[0] = i128::MAX -10000000;

    let mut info = -1;

    unsafe {
        ab13dd_(
            time_domain.as_ptr(),
            e_shape.as_ptr(),
            equilibration.as_ptr(),
            d_is_nonzero.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            f_peak.as_mut_ptr(),
            a.as_mut_ptr(),
            &(a.nrows() as c_int),
            e.as_mut_ptr(),
            &(e.nrows() as c_int),
            b.as_mut_ptr(),
            &(b.nrows() as c_int),
            c.as_mut_ptr(),
            &(c.nrows() as c_int),
            d.as_mut_ptr(),
            &(d.nrows() as c_int),
            peak_gain.as_mut_ptr(),
            &tol,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            cwork.as_mut_ptr(),
            &(lcwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(format!(
            "SLICOT ab13dd returned with error code. info = {}",
            info
        ));
    }

    let norm = peak_gain[0];
    Ok(norm)
}

/// Computes the invariant zeros of a continuous-time state-space system.
///
/// # Arguments
/// * `a` - State matrix \( A \)
/// * `b` - Input matrix \( B \)
/// * `c` - Output matrix \( C \)
/// * `d` - Feedthrough matrix \( D \)
///
/// # Returns
/// A `Result` containing a vector of complex numbers representing the system
/// zeros, or an error message.
///
/// Uses SLICOT and LAPACK for numerical computations.
pub fn zeros(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
) -> Result<Vec<Complex64>, String> {
    Ss::<Continuous>::verify_dimensions(a, b, c, d)
        .map_err(|e| e.to_string())?;

    let scaling = CString::new("S").unwrap();

    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();

    let mut num_inv_zeros = -1 as c_int;
    let mut rank_tf = -1 as c_int;
    let mut degree_inf_div = -1 as c_int;
    let mut num_right_kronecker_indecies = -1 as c_int;
    let mut num_left_kronecker_indecies = -1 as c_int;

    let mut inf_elementrary_divisors = vec![-1 as c_int; n];
    let mut right_kronecker_indecies = vec![-1 as c_int; n.max(m) + 1];
    let mut left_kronecker_indecies = vec![-1 as c_int; n.max(p) + 1];

    let mut a_f = DMatrix::zeros(n + m, n + p.min(m));
    let mut b_f = DMatrix::zeros(n + p, n + m);

    let tol = 1e-3 as c_double;

    let mut iwork = vec![0 as c_int; m.max(p)];
    let mut dwork = vec![0.0 as c_double; 1];
    let ldwork = -1 as c_int; // query

    let mut info = -1 as c_int;
    // query workspace size
    unsafe {
        ab08nd_(
            scaling.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            a.as_ptr(),
            &(a.nrows() as c_int),
            b.as_ptr(),
            &(b.nrows() as c_int),
            c.as_ptr(),
            &(c.nrows() as c_int),
            d.as_ptr(),
            &(d.nrows() as c_int),
            &mut num_inv_zeros,
            &mut rank_tf,
            &mut degree_inf_div,
            &mut num_right_kronecker_indecies,
            &mut num_left_kronecker_indecies,
            inf_elementrary_divisors.as_mut_ptr(),
            right_kronecker_indecies.as_mut_ptr(),
            left_kronecker_indecies.as_mut_ptr(),
            a_f.as_mut_ptr(),
            &(a_f.nrows() as c_int),
            b_f.as_mut_ptr(),
            &(b_f.nrows() as c_int),
            &tol,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(format!(
            "SLICOT failed to find zeros of state-space system. Info = {}",
            info
        ));
    }
    let ldwork = dwork[0] as usize;
    assert!(ldwork > 0);
    dwork = vec![0.0 as c_double; ldwork];
    unsafe {
        ab08nd_(
            scaling.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            a.as_ptr(),
            &(a.nrows() as c_int),
            b.as_ptr(),
            &(b.nrows() as c_int),
            c.as_ptr(),
            &(c.nrows() as c_int),
            d.as_ptr(),
            &(d.nrows() as c_int),
            &mut num_inv_zeros,
            &mut rank_tf,
            &mut degree_inf_div,
            &mut num_right_kronecker_indecies,
            &mut num_left_kronecker_indecies,
            inf_elementrary_divisors.as_mut_ptr(),
            right_kronecker_indecies.as_mut_ptr(),
            left_kronecker_indecies.as_mut_ptr(),
            a_f.as_mut_ptr(),
            &(a_f.nrows() as c_int),
            b_f.as_mut_ptr(),
            &(b_f.nrows() as c_int),
            &tol,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(format!(
            "SLICOT failed to find zeros of state-space system. Info = {}",
            info
        ));
    }
    if num_inv_zeros == 0 {
        return Ok(vec![]);
    }

    /////////////////////////////////
    // LAPACK generalized eigenvalue problem
    ////////////////////////////////////
    let comp_left_eigen_vectors = CString::new("N").unwrap();
    let comp_right_eigen_vectors = CString::new("N").unwrap();

    let n_gen_eigen = num_inv_zeros as usize;
    let mut alpha_re = vec![0.0 as c_double; n_gen_eigen];
    let mut alpha_im = vec![0.0 as c_double; n_gen_eigen];
    let mut beta = vec![0.0 as c_double; n_gen_eigen];

    let mut eigen_vec_left = DMatrix::zeros(n_gen_eigen, n_gen_eigen);
    let mut eigen_vec_right = DMatrix::zeros(n_gen_eigen, n_gen_eigen);
    let mut work = vec![0.0 as c_double; 1];
    let mut lwork = -1 as c_int;
    let mut info = -1 as c_int;

    // query workspace size
    unsafe {
        dggev_(
            comp_left_eigen_vectors.as_ptr(),
            comp_right_eigen_vectors.as_ptr(),
            &(n_gen_eigen as c_int),
            a_f.as_mut_ptr(),
            &(a_f.nrows() as c_int),
            b_f.as_mut_ptr(),
            &(b_f.nrows() as c_int),
            alpha_re.as_mut_ptr(),
            alpha_im.as_mut_ptr(),
            beta.as_mut_ptr(),
            eigen_vec_left.as_mut_ptr(),
            &(n_gen_eigen as c_int),
            eigen_vec_right.as_mut_ptr(),
            &(n_gen_eigen as c_int),
            work.as_mut_ptr(),
            &lwork,
            &mut info,
        );
    }
    if info != 0 {
        return Err(format!("LAPACK DGGEV returned info = {}", info));
    }
    lwork = dwork[0] as c_int;
    assert!(lwork > 0);
    work = vec![0.0; lwork as usize];

    // Calculate
    unsafe {
        dggev_(
            comp_left_eigen_vectors.as_ptr(),
            comp_right_eigen_vectors.as_ptr(),
            &(n_gen_eigen as c_int),
            a_f.as_mut_ptr(),
            &(a_f.nrows() as c_int),
            b_f.as_mut_ptr(),
            &(b_f.nrows() as c_int),
            alpha_re.as_mut_ptr(),
            alpha_im.as_mut_ptr(),
            beta.as_mut_ptr(),
            eigen_vec_left.as_mut_ptr(),
            &(n_gen_eigen as c_int),
            eigen_vec_right.as_mut_ptr(),
            &(n_gen_eigen as c_int),
            work.as_mut_ptr(),
            &lwork,
            &mut info,
        );
    }
    if info != 0 {
        return Err(format!("LAPACK DGGEV returned info = {}", info));
    }

    let mut zeros = Vec::with_capacity(n_gen_eigen);
    for i in 0..n_gen_eigen {
        let den = beta[i];
        assert!(den != 0.0);
        let zero = c64(alpha_re[i] / den, alpha_im[i] / den);
        zeros.push(zero);
    }
    Ok(zeros)
}

////////////////////////////////////////////////////////////
// State-Space methods
////////////////////////////////////////////////////////////

impl<U: Time + 'static> Ss<U> {
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
    pub fn norm_h2(&self) -> Result<f64, String> {
        h2_norm::<U>(
            self.a().clone(),
            self.b().clone(),
            self.c().clone(),
            self.d().clone(),
        )
        .map_err(|e| e.to_string())
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
    pub fn norm_hinf(&self) -> Result<f64, String> {
        hinf_norm::<U>(
            self.a().clone(),
            self.b().clone(),
            self.c().clone(),
            self.d().clone(),
        )
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

    /// Checks if the system is stable
    ///
    /// The system is stable if all poles are in the left half plane. I.e.
    /// Re(pole) < 0.
    ///
    /// # Returns
    /// `bool` true if the system is stable.
    pub fn is_stable(&self) -> bool {
        let poles = self.poles();
        poles.iter().all(|pole| pole.re < 0.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        systems::Tf,
        transformations::SsRealization::{ControllableCF, ObservableCF},
    };

    use approx::assert_abs_diff_eq;
    use num_complex::c64;
    use rand::Rng;
    #[test]
    fn ss_system_norms() {
        let sys_tf = Tf::s() / (Tf::s() + 1.0);
        let sys = sys_tf.to_ss_method(ObservableCF).unwrap();
        assert!(sys.clone().norm_h2().is_err());
        assert_abs_diff_eq!(sys.clone().norm_hinf().unwrap(), 1.0);

        let sys = 2.0 / (Tf::s() + 1.0);
        let sys = sys.to_ss_method(ObservableCF).unwrap();
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
            let sys = sys_tf.to_ss_method(ObservableCF).unwrap();
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
                    sys_tf.to_ss_method(ObservableCF).unwrap().poles();

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

        let ss = (1.0 / Tf::s()).to_ss_method(ObservableCF).unwrap();
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
        let zeros = tf.to_ss_method(ControllableCF).unwrap().zeros().unwrap();
        println!("zeros: {:?}", zeros);
        assert_eq!(zeros.len(), 2);

        let tf = 1.0 / Tf::s();
        let zeros = tf.to_ss_method(ObservableCF).unwrap().zeros().unwrap();
        assert_eq!(zeros.len(), 0);

        for _ in 0..100 {
            let mut rng = rand::rng();
            for zero in (0..100).map(|_| rng.random_range(-100.0..100.0)) {
                let tf = (Tf::s() - zero) / (Tf::s() + 1.0).powi(2);
                let sys = tf.to_ss_method(ObservableCF).unwrap();
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
                let sys = tf.to_ss_method(ObservableCF).unwrap();
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

    #[test]
    fn ss_is_stable() {
        let tf = 1.0 / Tf::s();
        assert_eq!(tf.to_ss_method(ObservableCF).unwrap().is_stable(), false);

        let tf = (Tf::s() + 1.0) / (Tf::s() + 1.0).powi(4);
        let ss = tf.to_ss_method(ObservableCF).unwrap();
        assert_eq!(ss.is_stable(), true);

        let tf = 1.0 / ((Tf::s() + 1.0) * (Tf::s() - 2.0));
        let ss = tf.to_ss_method(ObservableCF).unwrap();
        assert_eq!(ss.is_stable(), false);

        let tf = 1.0 / (Tf::s().powi(2) + 2.0 * 0.01 * Tf::s() + 1.0);
        let ss = tf.to_ss_method(ObservableCF).unwrap();
        assert_eq!(ss.is_stable(), true);
    }
}
