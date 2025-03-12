extern crate libc;
use libc::{c_char, c_double, c_long};

use nalgebra::DMatrix;

use crate::{ss::Ss, tf::Tf, traits::Time};
use std::{error::Error, ffi::CString};

#[link(name = "slicot")]
#[link(name = "blas")]
#[link(name = "lapack")]
unsafe extern "C" {
    fn tb04ad_(
        rowcol: *const c_char, // CHARACTER*1 -> Pointer to a single char
        n: *const c_long,      // INTEGER
        m: *const c_long,
        p: *const c_long,
        a: *mut c_double, // DOUBLE PRECISION array
        lda: *const c_long,
        b: *mut c_double,
        ldb: *const c_long,
        c: *mut c_double,
        ldc: *const c_long,
        d: *const c_double,
        ldd: *const c_long,
        nr: *mut c_long,       // Output INTEGER
        index: *mut c_long,    // INTEGER array
        dcoeff: *mut c_double, // DOUBLE PRECISION array
        lddcoe: *const c_long,
        ucoeff: *mut c_double, // DOUBLE PRECISION array
        lduco1: *const c_long,
        lduco2: *const c_long,
        tol1: *const c_double,
        tol2: *const c_double,
        iwork: *mut c_long,   // INTEGER array
        dwork: *mut c_double, // DOUBLE PRECISION array
        ldwork: *const c_long,
        info: *mut c_long, // Output INTEGER
    );
}

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
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("SLICOT tb04ad_ failed with info code {}", info),
        )));
    }

    let den_degree = index[0] as usize;
    let den: Vec<f64> = dcoeff.view((0, 0), (1, den_degree + 1)).iter().rev().map(|&x| x).collect();
    let num: Vec<f64> = ucoeff.view((0, 0), (1, den_degree + 1)).iter().rev().map(|&x| x).collect();
    let tf = Tf::<f64, U>::new(num.as_slice(), den.as_slice());
    println!("tf:\n{}", tf);

    Ok(tf)
}
