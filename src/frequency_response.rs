use num_complex::{Complex64, c64};

use crate::{
    tf::Tf,
    traits::{Mag2Db, Rad2Deg, Time},
};

pub fn lin_space(
    start: f64,
    end: f64,
    n: usize,
) -> impl ExactSizeIterator<Item = f64> {
    assert!(n > 1, "n must be greater than one");
    let step = (end - start) / (n as f64 - 1.0);
    (0..n).map(move |i| start + step * i as f64)
}

pub fn log_space(
    start: f64,
    end: f64,
    n: usize,
    base: usize,
) -> impl ExactSizeIterator<Item = f64> {
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

    nums.map(move |x| (base as f64).powf(x))
}

pub fn bode<U: Time>(
    sys: Tf<f64, U>,
    min_freq: f64,
    max_freq: f64,
) -> Vec<[f64; 3]> {
    let freqs = log_space(min_freq, max_freq, 1000, 10);
    bode_freqs(sys, freqs)
}

pub fn bode_freqs<U: Time>(
    sys: Tf<f64, U>,
    freqs: impl Iterator<Item = f64>,
) -> Vec<[f64; 3]> {
    let mut mag_phase_freq_vec = Vec::with_capacity(freqs.size_hint().0);

    for omega in freqs {
        let c = c64(0., omega);
        let sys_val = sys.eval(&c);
        mag_phase_freq_vec.push([
            sys_val.norm().mag2db(),
            sys_val.arg().rad2deg(),
            omega,
        ]);
    }
    mag_phase_freq_vec
}

pub fn nyquist<U: Time>(
    sys: Tf<f64, U>,
    min_freq: f64,
    max_freq: f64,
) -> Vec<Complex64> {
    let freqs = log_space(min_freq, max_freq, 1000, 10);
    nyquist_freqs(sys, freqs)
}

pub fn nyquist_freqs<U: Time>(
    sys: Tf<f64, U>,
    freqs: impl Iterator<Item = f64>,
) -> Vec<Complex64> {
    let mut pos_vals = Vec::with_capacity(freqs.size_hint().0);
    let mut neg_vals = Vec::with_capacity(freqs.size_hint().0);

    for freq in freqs {
        pos_vals.push(sys.eval(&c64(0., freq)));
        neg_vals.push(sys.eval(&c64(0., -freq)));
    }

    pos_vals.extend(neg_vals.iter().rev());
    pos_vals
}
