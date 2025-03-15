use num_complex::{Complex64, c64};

use crate::{
    tf::Tf,
    traits::{Mag2Db, Rad2Deg, Time},
};

pub fn lin_space(start: f64, end: f64, n: usize) -> Vec<f64> {
    assert!(n > 1, "n must be greater than one");
    let step = (end - start) / (n as f64 - 1.0);
    (0..n).map(|i| start + step * i as f64).collect()
}

pub fn log_space(start: f64, end: f64, n: usize, base: usize) -> Vec<f64> {
    assert!(
        start > 0.,
        "logarithm of negative numbers are not implemented"
    );
    assert!(
        end > 0.,
        "logarithm of negative numbers are not implemented"
    );
    let start_log = start.log(base as f64);
    let end_log = end.log(base as f64);

    let nums = lin_space(start_log, end_log, n);

    nums.iter().map(|x| (base as f64).powf(*x)).collect()
}

pub fn bode<U: Time>(
    sys: Tf<f64, U>,
    min_freq: f64,
    max_freq: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let freqs = log_space(min_freq, max_freq, 1000, 10);
    let (mag_vec, phase_vec) = bode_freqs(sys, freqs.as_slice());
    (mag_vec, phase_vec, freqs)
}

pub fn bode_freqs<U: Time>(
    sys: Tf<f64, U>,
    freqs: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let mut mag_vec = Vec::with_capacity(freqs.len());
    let mut phase_vec = Vec::with_capacity(freqs.len());

    for &omega in freqs {
        let c = c64(0., omega);
        let sys_val = sys.eval(&c);
        mag_vec.push(sys_val.norm().mag2db());
        phase_vec.push(sys_val.arg().rad2deg());
    }
    (mag_vec, phase_vec)
}

pub fn nyquist<U: Time>(
    sys: Tf<f64, U>,
    min_freq: f64,
    max_freq: f64,
) -> Vec<Complex64> {
    let freqs = log_space(min_freq, max_freq, 1000, 10);
    nyquist_freqs(sys, freqs.as_slice())
}

pub fn nyquist_freqs<U: Time>(
    sys: Tf<f64, U>,
    freqs: &[f64],
) -> Vec<Complex64> {
    let pos_freqs = freqs.iter().map(|freq| sys.eval(&c64(0., *freq)));
    let neg_freqs = freqs.iter().rev().map(|freq| sys.eval(&c64(0., -*freq)));

    pos_freqs.chain(neg_freqs).collect()
}
