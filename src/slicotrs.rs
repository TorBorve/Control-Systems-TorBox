extern crate libc;
extern crate netlib_src; // enable linking with blas and lapack
use crate::{Continuous, Discrete, ss::Ss, tf::Tf, traits::Time};
use libc::{c_char, c_double, c_int};
use nalgebra::DMatrix;
use std::{error::Error, ffi::CString};

unsafe extern "C" {
    /// SLICOT transform state-space to transfer matrix. See [SLICOT documentation](https://github.com/TorBorve/SLICOT-Reference)
    fn tb04ad_(
        rowcol: *const c_char,
        n: *const c_int,
        m: *const c_int,
        p: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        c: *mut c_double,
        ldc: *const c_int,
        d: *const c_double,
        ldd: *const c_int,
        nr: *mut c_int,
        index: *mut c_int,
        dcoeff: *mut c_double,
        lddcoe: *const c_int,
        ucoeff: *mut c_double,
        lduco1: *const c_int,
        lduco2: *const c_int,
        tol1: *const c_double,
        tol2: *const c_double,
        iwork: *mut c_int,
        dwork: *mut c_double,
        ldwork: *const c_int,
        info: *mut c_int,
    );

    // Calculate H2 or L2 norm of State space system
    fn ab13bd_(
        dico: *const c_char,
        jobn: *const c_char,
        n: *const c_int,
        m: *const c_int,
        p: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        c: *mut c_double,
        ldc: *const c_int,
        d: *mut c_double,
        ldd: *const c_int,
        nq: *mut c_int,
        tol: *const c_double,
        dwork: *mut c_double,
        ldwork: *const c_int,
        iwarn: *mut c_int,
        info: *mut c_int,
    ) -> c_double;

    // Calculate H-inf norm of State space system
    fn ab13dd_(
        dico: *const c_char,
        jobe: *const c_char,
        equil: *const c_char,
        jobd: *const c_char,
        n: *const c_int,
        m: *const c_int,
        p: *const c_int,
        fpeak: *mut c_double,
        a: *mut c_double,
        lda: *const c_int,
        e: *mut c_double,
        lde: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        c: *mut c_double,
        ldc: *const c_int,
        d: *mut c_double,
        ldd: *const c_int,
        gpeak: *mut c_double,
        tol: *const c_double,
        iwork: *mut c_int,
        dwork: *mut c_double,
        ldwork: *const c_int,
        cwork: *mut c_double,
        lcwork: *const c_int,
        info: *mut c_int,
    );
}

/// Converts a state-space system (A, B, C, D) into a transfer function using
/// the SLICOT `tb04ad_` function.
///
/// # Parameters:
/// - `a`: A `DMatrix` representing the state matrix `A`.
/// - `b`: A `DMatrix` representing the input matrix `B`.
/// - `c`: A `DMatrix` representing the output matrix `C`.
/// - `d`: A `DMatrix` representing the direct transmission matrix `D`.
///
/// # Returns:
/// A `Result` containing the transfer function (`Tf<f64, U>`) or an error
/// (`Box<dyn Error>`).
pub fn ss2tf_tb04ad<U: Time + 'static>(
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
pub fn h2_norm<U: Time + 'static>(
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
pub fn hinf_norm<U: Time + 'static>(
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
