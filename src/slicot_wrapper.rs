use std::ffi::{c_char, c_double, c_int};

unsafe extern "C" {
    /// SLICOT transform state-space to transfer matrix. See [SLICOT documentation](https://github.com/TorBorve/SLICOT-Reference)
    pub fn tb04ad_(
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
}

unsafe extern "C" {

    // Calculate H2 or L2 norm of State space system
    pub fn ab13bd_(
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
    pub fn ab13dd_(
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

    // Create regular pencil for system, which can be used to find invariant
    // zeros
    pub fn ab08nd_(
        equil: *const c_char,
        n: *const c_int,
        m: *const c_int,
        p: *const c_int,
        a: *const c_double,
        lda: *const c_int,
        b: *const c_double,
        ldb: *const c_int,
        c: *const c_double,
        ldc: *const c_int,
        d: *const c_double,
        ldd: *const c_int,
        nu: *mut c_int,
        rank: *mut c_int,
        dinfz: *mut c_int,
        nkror: *mut c_int,
        nkrol: *mut c_int,
        infz: *mut c_int,
        kronr: *mut c_int,
        kronl: *mut c_int,
        a_f: *mut c_double,
        ldaf: *const c_int,
        b_f: *mut c_double,
        ldbf: *const c_int,
        tol: *const c_double,
        iwork: *mut c_int,
        dwork: *mut c_double,
        ldwork: *const c_int,
        info: *mut c_int,
    );
}

// LAPACK function
unsafe extern "C" {
    // Compute eigenvalues and vectors for generalized eigenvalue problem
    pub fn dggev_(
        jobvl: *const c_char,
        jobvr: *const c_char,
        n: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        alphar: *mut c_double,
        alphai: *mut c_double,
        beta: *mut c_double,
        vl: *mut c_double,
        ldvl: *const c_int,
        vr: *mut c_double,
        ldvr: *const c_int,
        work: *mut c_double,
        lwork: *const c_int,
        info: *mut c_int,
    );
}
