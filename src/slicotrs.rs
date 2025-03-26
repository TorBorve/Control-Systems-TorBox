extern crate libc;
extern crate netlib_src; // enable linking with blas and lapack
use crate::{Continuous, ss::Ss, tf::Tf, traits::Time};
use libc::{c_char, c_double, c_long};
use nalgebra::DMatrix;
use std::{error::Error, ffi::CString};

// #[derive(Clone, Copy, Debug)]
// #[repr(C)]
// struct ComplexMock {
//     pub re: f64,
//     pub im: f64,
// }

unsafe extern "C" {
    /// SLICOT transform state-space to transfer matrix. See [SLICOT documentation](https://github.com/TorBorve/SLICOT-Reference)
    fn tb04ad_(
        rowcol: *const c_char,
        n: *const c_long,
        m: *const c_long,
        p: *const c_long,
        a: *mut c_double,
        lda: *const c_long,
        b: *mut c_double,
        ldb: *const c_long,
        c: *mut c_double,
        ldc: *const c_long,
        d: *const c_double,
        ldd: *const c_long,
        nr: *mut c_long,
        index: *mut c_long,
        dcoeff: *mut c_double,
        lddcoe: *const c_long,
        ucoeff: *mut c_double,
        lduco1: *const c_long,
        lduco2: *const c_long,
        tol1: *const c_double,
        tol2: *const c_double,
        iwork: *mut c_long,
        dwork: *mut c_double,
        ldwork: *const c_long,
        info: *mut c_long,
    );

    // Calculate H2 or L2 norm of State space system
    fn ab13bd_(
        dico: *const c_char,
        jobn: *const c_char,
        n: *const c_long,
        m: *const c_long,
        p: *const c_long,
        a: *mut c_double,
        lda: *const c_long,
        b: *mut c_double,
        ldb: *const c_long,
        c: *mut c_double,
        ldc: *const c_long,
        d: *mut c_double,
        ldd: *const c_long,
        nq: *mut c_long,
        tol: *const c_double,
        dwork: *mut c_double,
        ldwork: *const c_long,
        iwarn: *mut c_long,
        info: *mut c_long,
    ) -> c_double;

    // Calculate H-inf norm of State space system
    fn ab13dd_(
        dico: *const c_char,
        jobe: *const c_char,
        equil: *const c_char,
        jobd: *const c_char,
        n: *const c_long,
        m: *const c_long,
        p: *const c_long,
        fpeak: *mut c_double,
        a: *mut c_double,
        lda: *const c_long,
        e: *mut c_double,
        lde: *const c_long,
        b: *mut c_double,
        ldb: *const c_long,
        c: *mut c_double,
        ldc: *const c_long,
        d: *mut c_double,
        ldd: *const c_long,
        gpeak: *mut c_double,
        tol: *const c_double,
        iwork: *mut c_long,
        dwork: *mut c_double,
        ldwork: *const c_long,
        cwork: *mut u128,
        lcwork: *const c_long,
        info: *mut c_long,
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

    let mut nr: c_long = 0;
    let mut index = vec![0 as c_long; p];

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
    let mut iwork = vec![0 as c_long; n + max(m, p)];

    let mp = m;
    let pm = p;
    let ldwork = 10
        * max(
            1,
            n * (n + 1) + max(max(n * mp + 2 * n + max(n, mp), 3 * mp), pm),
        );
    let mut dwork = vec![0. as c_double; ldwork];
    let mut info = 0 as c_long;

    unsafe {
        tb04ad_(
            rowcol.as_ptr(),
            &(n as c_long),
            &(m as c_long),
            &(p as c_long),
            a_in.as_mut_ptr(),
            &(lda as c_long),
            b_in.as_mut_ptr(),
            &(ldb as c_long),
            c_in.as_mut_ptr(),
            &(ldc as c_long),
            d_in.as_mut_ptr(),
            &(ldd as c_long),
            &mut nr,
            index.as_mut_ptr(),
            dcoeff.as_mut_ptr(),
            &(lddcoe as c_long),
            ucoeff.as_mut_ptr(),
            &(lduco1 as c_long),
            &(lduco2 as c_long),
            &tol1,
            &tol2,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_long),
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
    } else if std::any::TypeId::of::<U>()
        == std::any::TypeId::of::<Continuous>()
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

    let mut nq = -1 as c_long;
    let tol = -1 as c_double;

    let ldwork = 5 * 1
        .max(m * (n + m) + (n * (n + 5)).max((m * (m + 2)).max(4 * p)))
        .max(n * (n.max(p) + 4) + n.min(p));
    let mut dwork = vec![0.0 as c_double; ldwork];

    let mut iwarn = -1 as c_long;
    let mut info = -1 as c_long;

    let h2_norm = unsafe {
        ab13bd_(
            time_domain.as_ptr(),
            h2_or_l2.as_ptr(),
            &(n as c_long),
            &(m as c_long),
            &(p as c_long),
            a.as_mut_ptr(),
            &(lda as c_long),
            b.as_mut_ptr(),
            &(ldb as c_long),
            c.as_mut_ptr(),
            &(ldc as c_long),
            d.as_mut_ptr(),
            &(ldd as c_long),
            &mut nq,
            &tol,
            dwork.as_mut_ptr(),
            &(ldwork as c_long),
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
    } else if std::any::TypeId::of::<U>()
        == std::any::TypeId::of::<Continuous>()
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

    let mut f_peak = [0.0, 1.0]; // initial guess not active

    let mut e = DMatrix::identity(n, n);

    let mut peak_gain = [-1.0, -1.0];
    let tol = 0.0;

    let mut iwork = vec![0; 2*n]; // Wihtout multiplication with 2 the code fails(?!)
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
    let mut dwork = vec![-1.0; ldwork];

    let lcwork = 10 * 1.max((n + m) * (n + p) + 2 * p.min(m) + p.max(m));
    let mut cwork = vec![1; 20 * (lcwork + 10)]; // Mulipy by two because cwork is complex 16 bytes not f64
    // cwork[0] = i128::MAX -10000000;

    let mut info = -1;

    unsafe {
        ab13dd_(
            time_domain.as_ptr(),
            e_shape.as_ptr(),
            equilibration.as_ptr(),
            d_is_nonzero.as_ptr(),
            &(n as c_long),
            &(m as c_long),
            &(p as c_long),
            f_peak.as_mut_ptr(),
            a.as_mut_ptr(),
            &(n as c_long),
            e.as_mut_ptr(),
            &(n as c_long),
            b.as_mut_ptr(),
            &(n as c_long),
            c.as_mut_ptr(),
            &(p as c_long),
            d.as_mut_ptr(),
            &(p as c_long),
            peak_gain.as_mut_ptr(),
            &tol,
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_long),
            cwork.as_mut_ptr().add(10), /* don't know why add is necessary. think there is some
                                        * issue with byte alignment or
                                        * similar... */
            &(lcwork as c_long),
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

#[cfg(test)]
mod tests {

    use crate::{Continuous, Tf, slicotrs::hinf_norm, tf2ss};

    use super::h2_norm;

    #[test]
    fn tets_h2_norm() {
        let sys_tf = 0.1 / (Tf::s() + 1.0);
        // let sys_tf = 1.0 / (Tf::s().powi(2) + 2.0 * 0.1 * Tf::s() + 2.0);
        let sys =
            tf2ss(sys_tf.clone(), crate::SsRealization::ObservableCF).unwrap();

        let h2norm = h2_norm::<Continuous>(
            sys.a().clone(),
            sys.b().clone(),
            sys.c().clone(),
            sys.d().clone(),
        )
        .unwrap();
        println!("H2 norm: {} of sys: \n{}", h2norm, sys_tf);

        let hinfnorm = hinf_norm::<Continuous>(
            sys.a().clone(),
            sys.b().clone(),
            sys.c().clone(),
            sys.d().clone(),
        )
        .unwrap();
        println!("H-inf norm: {}", hinfnorm);
    }
}
