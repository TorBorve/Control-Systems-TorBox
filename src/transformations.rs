use nalgebra::DMatrix;
use std::error::Error;

use crate::{
    Continuous,
    slicot_wrapper::{tb01pd_, tb04ad_},
    systems::{Ss, Tf},
    utils::traits::Time,
};
use std::ffi::{CString, c_double, c_int};

/// Converts state-space matrices to a transfer function.
///
/// Uses SLICOT TB04AB
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
fn ss2tf_mat<U: Time + 'static>(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
) -> Result<Tf<f64, U>, Box<dyn Error + 'static>> {
    Ss::<U>::verify_dimensions(a, b, c, d)?;

    let rowcol = CString::new("R")?;

    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();

    if m != 1 || p != 1 {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "State Space system must be SISO single-input-single-output",
        )));
    }

    let lda = n;
    let ldb = n;
    let ldc = p;
    let ldd = p;

    let mut a_in = DMatrix::zeros(lda, n);
    a_in.view_mut((0, 0), (n, n)).copy_from(a);
    let mut b_in = DMatrix::zeros(ldb, m);
    b_in.view_mut((0, 0), (n, m)).copy_from(b);
    let mut c_in = DMatrix::zeros(ldc, n);
    c_in.view_mut((0, 0), (p, n)).copy_from(c);
    let mut d_in = DMatrix::zeros(ldd, m);
    d_in.view_mut((0, 0), (p, m)).copy_from(d);

    let mut nr: c_int = 0;
    let mut index = vec![0 as c_int; p];

    let lddcoe = p;
    let mut dcoeff = DMatrix::<f64>::zeros(lddcoe, n + 1);

    let lduco1 = p;
    let lduco2 = m;
    assert!(p == 1 && m == 1);
    let mut ucoeff = DMatrix::<f64>::zeros(lduco1, n + 1);

    // use default tolerances
    let tol1 = -1. as c_double;
    let tol2 = -1 as c_double;

    use std::cmp::max;
    let mut iwork = vec![0 as c_int; n + max(m, p)];

    let mp = m;
    let pm = p;
    let ldwork = 10
        * max(
            1,
            n * (n + 1) + max(max(n * mp + 2 * n + max(n, mp), 3 * mp), pm),
        );
    let mut dwork = vec![0. as c_double; ldwork];
    let mut info = 0 as c_int;

    unsafe {
        tb04ad_(
            rowcol.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            a_in.as_mut_ptr(),
            &(lda as c_int),
            b_in.as_mut_ptr(),
            &(ldb as c_int),
            c_in.as_mut_ptr(),
            &(ldc as c_int),
            d_in.as_mut_ptr(),
            &(ldd as c_int),
            &mut nr,
            index.as_mut_ptr(),
            dcoeff.as_mut_ptr(),
            &(lddcoe as c_int),
            ucoeff.as_mut_ptr(),
            &(lduco1 as c_int),
            &(lduco2 as c_int),
            &tol1,
            &tol2,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(Box::new(std::io::Error::other(format!(
            "SLICOT tb04ad_ failed with info code {}",
            info
        ))));
    }

    let den_degree = index[0] as usize;
    let den: Vec<f64> = dcoeff
        .view((0, 0), (1, den_degree + 1))
        .iter()
        .rev()
        .copied()
        .collect();
    let num: Vec<f64> = ucoeff
        .view((0, 0), (1, den_degree + 1))
        .iter()
        .rev()
        .copied()
        .collect();
    let tf = Tf::<f64, U>::new(num.as_slice(), den.as_slice());
    Ok(tf)
}

impl<U: Time + 'static> Tf<f64, U> {
    fn to_ss_observable(&self) -> Result<Ss<U>, Box<dyn Error + 'static>> {
        if !self.is_proper() {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Transfer function must be proper",
            )));
        }

        let tf = self.normalize();

        let (_, den_deg) = tf.degree_num_den();
        let nx = den_deg;

        let mut d = DMatrix::zeros(1, 1);
        if !tf.is_strictly_proper() {
            assert!(tf.numerator().coeffs().len() > nx);
            d[(0, 0)] = tf.numerator().coeffs()[nx];
        }

        if nx == 0 {
            assert_eq!(tf.relative_degree(), 0);
            return Ok(Ss::new_from_scalar(d[(0, 0)]));
        }

        let mut a = DMatrix::zeros(nx, nx);
        a.view_mut((1, 0), (nx - 1, nx - 1))
            .copy_from(&DMatrix::identity(nx - 1, nx - 1));
        for row in 0..nx {
            a[(row, nx - 1)] = -tf.denominator().coeffs()[row];
        }


        let mut num_extended = tf.numerator().coeffs().to_vec();
        num_extended.resize(nx, 0.);

        let mut b_values = num_extended;
        assert!(tf.denominator().coeffs().len() - 1 <= b_values.len());
        for (i, den_i) in tf.denominator().coeffs().iter().enumerate().take(nx)
        {
            b_values[i] += -d[(0, 0)] * den_i;
        }

        let b = DMatrix::from_column_slice(nx, 1, &b_values);

        let mut c = DMatrix::zeros(1, nx);
        c[(0, nx - 1)] = 1.;

        Ss::<U>::new(a, b, c, d)
    }

    fn to_ss_controllable(&self) -> Result<Ss<U>, Box<dyn Error + 'static>> {
        let ss = self.to_ss_observable()?;

        let a = ss.a().transpose();
        let b = ss.c().transpose();
        let c = ss.b().transpose();
        let d = ss.d().transpose();

        Ss::<U>::new(a, b, c, d)
    }

    /// Converts a transfer function to a state-space representation. Using
    /// Observable Canonical Form.
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing a state-space system (`Ss<U>`) or an error
    /// if the conversion fails.
    ///
    /// # Example
    ///
    /// ```rust
    /// use control_systems_torbox::*;
    /// let tf = 1.0/Tf::s();
    /// let ss = tf.to_ss();
    /// ```
    pub fn to_ss(&self) -> Result<Ss<U>, String> {
        self.to_ss_method(SsRealization::ObservableCF)
    }
    /// Converts a transfer function to a state-space representation.
    ///
    /// # Arguments
    ///
    /// * `method` - The realization method used for the conversion such as
    ///   Controllable Canonical Form or Observable Canonical Form.
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing a state-space system (`Ss<U>`) or an error
    /// if the conversion fails.
    ///
    /// # Example
    ///
    /// ```rust
    /// use control_systems_torbox::*;
    /// let tf = 1.0/Tf::s();
    /// let ss = tf.to_ss_method(SsRealization::ControllableCF);
    /// ```
    pub fn to_ss_method(&self, method: SsRealization) -> Result<Ss<U>, String> {
        match method {
            SsRealization::ObservableCF => {
                self.to_ss_observable().map_err(|e| e.to_string())
            }
            SsRealization::ControllableCF => {
                self.to_ss_controllable().map_err(|e| e.to_string())
            }
        }
    }
}

impl<U: Time + 'static> Ss<U> {
    /// Converts a state-space representation to a transfer function.
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing a transfer function (`Tf<f64, U>`) or an
    /// error if the conversion fails.
    ///
    /// # Example
    ///
    /// ```rust
    /// use control_systems_torbox::*;
    /// let ss = (1.0/Tf::s()).to_ss().unwrap();
    /// let tf = ss.to_tf().unwrap();
    /// ```
    pub fn to_tf(&self) -> Result<Tf<f64, U>, String> {
        ss2tf_mat(self.a(), self.b(), self.c(), self.d())
            .map_err(|e| e.to_string())
    }
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

/// Options for the minimal realization (`minreal`) operation.
pub struct MinrealOptions {
    /// Numerical tolerance for rank decisions. If zero, the algorithm chooses
    /// an appropriate value internally.
    pub tolerance: f64,
    /// Whether to remove uncontrollable states from the system.
    pub remove_uncontrollable: bool,
    /// Whether to remove unobservable states from the system.
    pub remove_unobservable: bool,
}

impl Default for MinrealOptions {
    /// Provides default options:
    /// - `tolerance`: 0.0 (auto-determined)
    /// - `remove_uncontrollable`: true
    /// - `remove_unobservable`: true
    fn default() -> Self {
        Self {
            tolerance: 0.0,
            remove_uncontrollable: true,
            remove_unobservable: true,
        }
    }
}

/// Performs in-place minimal realization of a continuous-time linear system in
/// state-space form, using the SLICOT routine `TB01PD`.
///
/// This function reduces the order of the system by eliminating uncontrollable
/// and/or unobservable states.
///
/// # Arguments
///
/// * `a` - Mutable reference to the state matrix `A` (n×n).
/// * `b` - Mutable reference to the input matrix `B` (n×m).
/// * `c` - Mutable reference to the output matrix `C` (p×n).
/// * `options` - Optional settings for the reduction algorithm.
///
/// # Returns
///
/// * `Ok(new_order)` - The new system order after reduction.
/// * `Err(msg)` - If the SLICOT routine fails.
///
/// # Errors
///
/// Returns an error if the underlying SLICOT call fails or if illegal matrix
/// dimensions are provided.
///
/// # Panics
///
/// Panics if matrices have inconsistent dimensions or invalid options are
/// passed (e.g., negative tolerance).
pub fn minreal_mat_mut(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    options: Option<MinrealOptions>,
) -> Result<usize, String> {
    let options = options.unwrap_or_default();
    let tolerance = options.tolerance;
    let remove_uncontrollable = options.remove_uncontrollable;
    let remove_unobservable = options.remove_unobservable;
    assert!(tolerance >= 0.0);

    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();

    assert!(a.is_square());
    assert_eq!(b.nrows(), n);
    assert_eq!(c.ncols(), n);

    let removal_kind = if remove_uncontrollable && remove_unobservable {
        CString::new("M").unwrap()
    } else if remove_uncontrollable {
        CString::new("C").unwrap()
    } else if remove_unobservable {
        CString::new("O").unwrap()
    } else {
        // why did they call the function?
        return Ok(n);
    };

    let apply_scaling = CString::new("S").unwrap();

    let mut reduced_order = -1 as c_int;
    // let tol = 0.0;
    let mut iwork = vec![0; n + m.max(p)];

    let ldwork = 100 * (n + n.max(3 * m).max(3 * p));
    let mut dwork = vec![0.0; ldwork];
    let mut info = -1 as c_int;

    unsafe {
        tb01pd_(
            removal_kind.as_ptr(),
            apply_scaling.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            a.as_mut_ptr(),
            &(a.nrows() as c_int),
            b.as_mut_ptr(),
            &(b.nrows() as c_int),
            c.as_mut_ptr(),
            &(c.nrows() as c_int),
            &mut reduced_order,
            &tolerance,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(format!(
            "Minreal failed. SLICOT TB01PD returned info = {}",
            info
        ));
    }

    assert!(reduced_order >= 0);
    Ok(reduced_order as usize)
}

/// Performs minimal realization of a state-space model by removing
/// uncontrollable and/or unobservable states. Returns a new set of
/// reduced-order matrices.
///
/// This is a safe wrapper around `minreal_mut`, performing the reduction
/// non-destructively by cloning the input matrices.
///
/// # Arguments
///
/// * `a` - State matrix `A`.
/// * `b` - Input matrix `B`.
/// * `c` - Output matrix `C`.
/// * `options` - Optional settings for the reduction algorithm.
///
/// # Returns
///
/// * `Ok([A_new, B_new, C_new])` - Reduced system matrices.
/// * `Err(msg)` - If the SLICOT routine fails.
///
/// # Example
///
/// ```rust
/// use control_systems_torbox::transformations::minreal_mat;
/// use nalgebra::dmatrix;
/// let a = dmatrix![1.0, 0.0;
///                  0.0, -0.2];
/// let b = dmatrix![0.0; 1.0];
/// let c = dmatrix![0.0, 1.0];
/// let [a_new, b_new, c_new] = minreal_mat(a, b, c, None).unwrap();
/// assert_eq!(a_new.nrows(), 1);
/// ```
pub fn minreal_mat(
    mut a: DMatrix<f64>,
    mut b: DMatrix<f64>,
    mut c: DMatrix<f64>,
    options: Option<MinrealOptions>,
) -> Result<[DMatrix<f64>; 3], String> {
    let new_nx = minreal_mat_mut(&mut a, &mut b, &mut c, options)?;
    let a_new = a.view((0, 0), (new_nx, new_nx)).into_owned();
    let b_new = b.view((0, 0), (new_nx, b.ncols())).into_owned();
    let c_new = c.view((0, 0), (c.nrows(), new_nx)).into_owned();
    Ok([a_new, b_new, c_new])
}

impl<U: Time + 'static> Ss<U> {
    /// Returns a minimal realization of this state-space system, eliminating
    /// uncontrollable and/or unobservable states.
    ///
    /// This method preserves the original system and returns a new
    /// reduced-order system.
    ///
    /// # Arguments
    ///
    /// * `options` - Optional settings for the minimal realization algorithm.
    ///
    /// # Returns
    ///
    /// * `Ok(Ss)` - A new state-space system with reduced state dimension.
    /// * `Err(msg)` - If the realization fails (e.g., due to invalid dimensions
    ///   or internal errors).
    ///
    /// # Example
    ///
    /// ```rust
    /// use control_systems_torbox::{Tf, Ss};
    /// let tf = Tf::s() / Tf::s().powi(2);
    /// let system = tf.to_ss().unwrap();
    /// let reduced = system.minreal(None).unwrap();
    /// ```
    pub fn minreal(
        &self,
        options: Option<MinrealOptions>,
    ) -> Result<Self, String> {
        let [a_new, b_new, c_new] = minreal_mat(
            self.a().clone(),
            self.b().clone(),
            self.c().clone(),
            options,
        )?;
        Ss::new(a_new, b_new, c_new, self.d().clone())
            .map_err(|e| e.to_string())
    }
}

impl Tf<f64, Continuous> {
    /// Returns a Pade approximation of a pure time delay `e^{-sT}` as a
    /// transfer function.
    ///
    /// The Pade approximation is a rational function (ratio of two polynomials)
    /// that approximates the exponential delay term in the Laplace domain.
    ///
    /// # Arguments
    ///
    /// * `time_delay` - The delay time `T` (must be non-negative).
    /// * `order` - The order `n` of the approximation. Higher orders improve
    ///   accuracy, especially at higher frequencies, but increase computational
    ///   complexity.
    ///
    /// # Returns
    ///
    /// A `Tf<f64, Continuous>` representing the rational approximation of the
    /// delay.
    ///
    /// # Panics
    ///
    /// Panics if `time_delay` is negative.
    ///
    /// # Examples
    ///
    /// ```
    /// use control_systems_torbox::Tf;
    /// let delay = Tf::pade(0.5, 3);
    /// println!{"{}", delay};
    /// ```
    ///
    /// # References
    /// - [RATIONAL APPROXIMATION OF TIME DELAY](https://www2.humusoft.cz/www/papers/tcp09/035_hanta.pdf)
    pub fn pade(time_delay: f64, order: usize) -> Tf<f64, Continuous> {
        pade_approx(time_delay, order)
    }
}

fn pade_approx(time_delay: f64, order: usize) -> Tf<f64, Continuous> {
    assert!(time_delay >= 0.0);

    if order == 0 || time_delay == 0.0 {
        return Tf::new_from_scalar(1.0);
    }

    let mut num_coeffs = Vec::with_capacity(order + 1);
    let mut den_coeffs = Vec::with_capacity(order + 1);

    num_coeffs.push(1.0);
    den_coeffs.push(1.0);

    for i in 1..=order {
        // Not efficient calculation of coeffs, however likely fast enough.
        num_coeffs.push(time_delay.powi(i as i32) * pade_coeff_num(i, order));
        den_coeffs.push(time_delay.powi(i as i32) * pade_coeff_den(i, order));
    }
    Tf::new(&num_coeffs, &den_coeffs).normalize()
}

fn pade_coeff_den(i: usize, n: usize) -> f64 {
    assert!(i <= n);

    let num = factorial(2 * n - i) * factorial(n);
    let den = factorial(2 * n) * factorial(i) * factorial(n - i);
    num as f64 / den as f64
}

fn pade_coeff_num(i: usize, n: usize) -> f64 {
    assert!(i <= n);

    let mut abs_coeff = pade_coeff_den(i, n);
    if i % 2 != 0 {
        abs_coeff *= -1.;
    }
    abs_coeff
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use rayon::prelude::*;

    use crate::utils::traits::Continuous;
    use approx::assert_abs_diff_eq;
    use rand::{Rng, seq::IteratorRandom};

    fn rand_proper_tf<U: Rng>(
        rng: &mut U,
        max_order: usize,
    ) -> Tf<f64, Continuous> {
        let den_order = rng.random_range(1..=max_order) as usize;
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
    fn pade_approx_test() {
        let sys = Tf::pade(2.5, 0);
        assert_abs_diff_eq!(sys, Tf::new_from_scalar(1.0));
        let sys = Tf::pade(0.0, 3);
        assert_abs_diff_eq!(sys, Tf::new_from_scalar(1.0));

        let sys = Tf::pade(2.5, 1);
        let sys_matlab = Tf::new(&[0.8, -1.0], &[0.8, 1.0]);
        assert_abs_diff_eq!(sys, sys_matlab, epsilon = 0.1);

        let sys = Tf::pade(2.5, 4);
        let sys_matlab = Tf::new(
            &[43.01, -53.76, 28.8, -8.0, 1.0],
            &[43.01, 53.76, 28.8, 8.0, 1.0],
        );
        assert_abs_diff_eq!(sys, sys_matlab, epsilon = 0.1);
    }

    #[test]
    #[should_panic]
    fn pade_approx_negative() {
        let _sys = Tf::pade(-1.0, 0);
    }

    #[test]
    fn ss2tf_test() {
        let a = DMatrix::from_row_slice(2, 2, &[0., 1., 0., 0.]);
        let b = DMatrix::from_row_slice(2, 1, &[0., 1.]);
        let c = DMatrix::from_row_slice(1, 2, &[1., 0.]);
        let d = DMatrix::zeros(1, 1);

        let ss = Ss::<Continuous>::new(a, b, c, d).unwrap();
        let tf = ss.to_tf().unwrap();

        let tf_ans = 1. / Tf::s().powi(2);
        println!("ss2Tf: \n{}", tf);
        assert_abs_diff_eq!(tf, tf_ans);
    }

    #[test]
    fn tf2ss_test() {
        let tf = 1. / Tf::s().powi(2);

        let ss = tf.to_ss().unwrap();
        let tf_ret = ss.to_tf().unwrap();
        let a_ans = DMatrix::from_row_slice(2, 2, &[0., 0., 1., 0.]);
        let b_ans = DMatrix::from_row_slice(2, 1, &[1., 0.]);
        let c_ans = DMatrix::from_row_slice(1, 2, &[0., 1.]);
        let d_ans = DMatrix::from_row_slice(1, 1, &[0.]);
        assert_abs_diff_eq!(ss.a(), &a_ans);
        assert_abs_diff_eq!(ss.b(), &b_ans);
        assert_abs_diff_eq!(ss.c(), &c_ans);
        assert_abs_diff_eq!(ss.d(), &d_ans);

        assert_abs_diff_eq!(tf_ret.normalize(), tf.normalize());

        let tf = Tf::<_, Continuous>::new_from_scalar(1.0);
        let ss = tf.to_ss().unwrap();
        assert_eq!(ss.order(), 0);
        assert_eq!(ss.d()[(0, 0)], 1.0);

        assert_eq!(tf.is_strictly_proper(), false);
        assert_eq!(tf.is_proper(), true);

    }

    #[test]
    fn convert_between_tf_and_ss() {
        let start = Instant::now();
        let num_iter = 10000;
        (0..num_iter).into_par_iter().for_each(|_| {
            let mut rng = rand::rng();
            let tf = rand_proper_tf(&mut rng, 20);
            let methods =
                [SsRealization::ControllableCF, SsRealization::ObservableCF];
            let ss = tf
                .to_ss_method(*methods.iter().choose(&mut rng).unwrap())
                .unwrap();
            let tf_ret = ss.to_tf().unwrap().normalize();

            let tf = tf.normalize();
            assert_abs_diff_eq!(tf, tf_ret, epsilon = 1e-3);
        });
        println!(
            "Time ss2tf(tf2ss()) transforms avg: {:?}",
            start.elapsed() / num_iter
        );
    }

    #[test]
    fn test_minreal() {
        let tf = (Tf::s() + 1.0) / (Tf::s().powi(2) + 2.0);
        let ss = tf.to_ss().unwrap();
        let tf_minreal = ss.minreal(None).unwrap().to_tf().unwrap();
        assert_abs_diff_eq!(tf, tf_minreal, epsilon = 1e-2);

        let tf = (Tf::s() + 1.0) / ((Tf::s() + 1.0) * Tf::s());
        let ss = tf.to_ss().unwrap();
        println!("a: {}, b: {}, c: {}, d: {}", ss.a(), ss.b(), ss.c(), ss.d());
        let tf_minreal = ss.minreal(None).unwrap().to_tf().unwrap();
        println!("tf minreal:\n{}", tf_minreal);
        assert_abs_diff_eq!(1.0 / Tf::s(), tf_minreal, epsilon = 1e-2);

        let tf = Tf::s().powi(5) / Tf::s().powi(6);
        let ss = tf.to_ss().unwrap();
        let tf_minreal = ss.minreal(None).unwrap().to_tf().unwrap();
        assert_abs_diff_eq!(tf_minreal, 1.0 / Tf::s(), epsilon = 1e-2);

        let tf = Tf::s() / (Tf::s() + 1e-3);
        let ss = tf.to_ss().unwrap();
        let tf_minreal = ss.minreal(None).unwrap().to_tf().unwrap();
        assert_abs_diff_eq!(tf, tf_minreal, epsilon = 1e-6);

        // let mut opts = MinrealOptions::default();
        // opts.tolerance = 0.01;
        // let tf_minreal_high_tol =
        // ss.minreal(Some(opts)).unwrap().to_tf().unwrap();
        // println!("tf high tol: \n{}", tf_minreal_high_tol);

        // assert_abs_diff_eq!(Tf::new_from_scalar(1.0), tf_minreal_high_tol);
    }
}
