use std::ffi::{CString, c_int};

use nalgebra::{DMatrix, zero};
use num_complex::{Complex64, c64};

use crate::{
    Continuous, Ss,
    slicot_wrapper::tb05ad_,
    systems::Tf,
    utils::traits::{Mag2Db, Rad2Deg, Time},
};
/// Generates a linearly spaced iterator between `start` and `end`, inclusive.
///
/// # Arguments
/// - `start`: The starting value of the sequence.
/// - `end`: The ending value of the sequence.
/// - `n`: The number of points to generate (must be greater than 1).
///
/// # Returns
/// - An iterator producing `n` evenly spaced `f64` values from `start` to
///   `end`.
///
/// # Panics
/// - Panics if `n` is less than or equal to 1.
pub fn lin_space(
    start: f64,
    end: f64,
    n: usize,
) -> impl ExactSizeIterator<Item = f64> {
    assert!(n >= 1, "n must be greater than or equal to one");
    let step = (end - start) / (n as f64 - 1.0);
    (0..n).map(move |i| start + step * i as f64)
}

/// Generates a logarithmically spaced iterator between `start` and `end`, using
/// the specified logarithmic base.
///
/// # Arguments
/// - `start`: The starting value of the sequence (must be greater than 0).
/// - `end`: The ending value of the sequence (must be greater than 0).
/// - `n`: The number of points to generate.
/// - `base`: The logarithmic base to use for spacing.
///
/// # Returns
/// - An iterator producing `n` logarithmically spaced `f64` values from `start`
///   to `end`.
///
/// # Panics
/// - Panics if `start` or `end` is less than or equal to 0.
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
    assert!(base > 1, "log_base() must be well defined");
    let start_log = start.log(base as f64);
    let end_log = end.log(base as f64);

    let nums = lin_space(start_log, end_log, n);

    nums.map(move |x| (base as f64).powf(x))
}

pub trait FrequencyResponse {
    fn freq_response(&self, freq: &Complex64) -> Complex64;
    /// Computes the Bode plot (magnitude and phase) for a transfer function
    /// over a frequency range.
    ///
    /// # Arguments
    /// - `sys`: The transfer function of the system to evaluate.
    /// - `min_freq`: The minimum frequency for the plot.
    /// - `max_freq`: The maximum frequency for the plot.
    ///
    /// # Returns
    /// - A vector of `[magnitude (dB), phase (degrees), frequency]` tuples for
    ///   each evaluated frequency.
    fn bode(&self, min_freq: f64, max_freq: f64) -> Vec<[f64; 3]> {
        let freqs = log_space(min_freq, max_freq, 1000, 10);
        self.bode_freqs(freqs)
    }

    /// Computes the Bode plot (magnitude and phase) for a transfer function
    /// over a given set of frequencies.
    ///
    /// # Arguments
    /// - `sys`: The transfer function of the system to evaluate.
    /// - `freqs`: An iterator of frequencies to evaluate the system at.
    ///
    /// # Returns
    /// - A vector of `[magnitude (dB), phase (degrees), frequency]` tuples for
    ///   each evaluated frequency.
    fn bode_freqs(&self, freqs: impl Iterator<Item = f64>) -> Vec<[f64; 3]> {
        let mut mag_phase_freq_vec = Vec::with_capacity(freqs.size_hint().0);

        for omega in freqs {
            let c = c64(0., omega);
            let sys_val = self.freq_response(&c);
            mag_phase_freq_vec.push([
                sys_val.norm().mag2db(),
                sys_val.arg().rad2deg(),
                omega,
            ]);
        }
        mag_phase_freq_vec
    }

    /// Computes the Nyquist plot for a transfer function over a frequency
    /// range.
    ///
    /// # Arguments
    /// - `sys`: The transfer function of the system to evaluate.
    /// - `min_freq`: The minimum frequency for the plot.
    /// - `max_freq`: The maximum frequency for the plot.
    ///
    /// # Returns
    /// - A vector of complex numbers representing the Nyquist plot.
    fn nyquist(&self, min_freq: f64, max_freq: f64) -> Vec<Complex64> {
        let freqs = log_space(min_freq, max_freq, 1000, 10);
        self.nyquist_freqs(freqs)
    }

    /// Computes the Nyquist plot for a transfer function over a given set of
    /// frequencies.
    ///
    /// # Arguments
    /// - `sys`: The transfer function of the system to evaluate.
    /// - `freqs`: An iterator of frequencies to evaluate the system at.
    ///
    /// # Returns
    /// - A vector of complex numbers representing the Nyquist plot.
    fn nyquist_freqs(
        &self,
        freqs: impl Iterator<Item = f64>,
    ) -> Vec<Complex64> {
        let mut pos_vals = Vec::with_capacity(freqs.size_hint().0);
        let mut neg_vals = Vec::with_capacity(freqs.size_hint().0);

        for freq in freqs {
            pos_vals.push(self.freq_response(&c64(0., freq)));
            neg_vals.push(self.freq_response(&c64(0., -freq)));
        }

        pos_vals.extend(neg_vals.iter().rev());
        pos_vals
    }
}

impl<U: Time> FrequencyResponse for Tf<f64, U> {
    fn freq_response(&self, freq: &Complex64) -> Complex64 {
        self.eval(freq)
    }
}

impl<U: Time + 'static> FrequencyResponse for Ss<U> {
    fn freq_response(&self, freq: &Complex64) -> Complex64 {
        assert_eq!(self.ninputs(), 1);
        assert_eq!(self.noutputs(), 1);

        let freq_resp = freq_response_ss_mat(
            *freq,
            &mut self.a().clone(),
            &mut self.b().clone(),
            &mut self.c().clone(),
            &mut self.d().clone(),
        )
        .unwrap();
        assert_eq!(freq_resp.nrows(), 1);
        assert_eq!(freq_resp.ncols(), 1);

        freq_resp[(0, 0)]
    }
}

fn freq_response_ss_mat(
    freq: num_complex::Complex64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
) -> Result<DMatrix<Complex64>, String> {
    // TODO: Should be possible to call it without specifying the time domain.
    Ss::<Continuous>::verify_dimensions(a, b, c, d).unwrap();
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    assert_eq!(m, 1, "Function is only tested for SISO systems");
    assert_eq!(p, 1, "Function is only tested for SISO systems");

    let baleig = CString::new("A").unwrap(); // balance A and compute condition number
    let inita = CString::new("G").unwrap(); // general A matrix

    let freq_in = crate::slicot_wrapper::Complex64 {
        re: freq.re,
        im: freq.im,
    };

    let mut rank_condition = -1.0;
    let zero_complex = crate::slicot_wrapper::Complex64 { re: 0.0, im: 0.0 };
    let mut freq_response = vec![zero_complex; p * m];

    let mut eigen_val_re = vec![0.0; n];
    let mut eigen_val_im = vec![0.0; n];

    let mut h_inv_times_b = vec![zero_complex; n * m];
    let l_h_inv_times_b = n;

    let mut iwork = vec![0; n];
    let ldwork = 10 * 1.max(n + n.max(m - 1).max(p - 1)) + 500;
    let mut dwork = vec![0.0; ldwork];

    let lzwork = 10 * 1.max(n * n + 2 * n) + 100;
    let mut zwork = vec![zero_complex; lzwork];

    let mut info = -1;

    unsafe {
        tb05ad_(
            baleig.as_ptr(),
            inita.as_ptr(),
            &(n as c_int),
            &(m as c_int),
            &(p as c_int),
            &freq_in,
            a.as_mut_ptr(),
            &(a.nrows() as c_int),
            b.as_mut_ptr(),
            &(b.nrows() as c_int),
            c.as_mut_ptr(),
            &(c.nrows() as c_int),
            &mut rank_condition,
            freq_response.as_mut_ptr(),
            &(p as c_int),
            eigen_val_re.as_mut_ptr(),
            eigen_val_im.as_mut_ptr(),
            h_inv_times_b.as_mut_ptr(),
            &(l_h_inv_times_b as c_int),
            iwork.as_mut_ptr(),
            dwork.as_mut_ptr(),
            &(ldwork as c_int),
            zwork.as_mut_ptr(),
            &(lzwork as c_int),
            &mut info,
        );
    }

    if info != 0 {
        return Err(format!(
            "Failed to compute frequency response of state space system. Slicot TB05AD failed with error code: {}",
            info
        ));
    }

    let mut response_matrix: DMatrix<Complex64> = DMatrix::zeros(p, m);

    for row in 0..p {
        for col in 0..m {
            let resp_no_d = freq_response[row + p * col];
            response_matrix[(row, col)] =
                c64(d[(row, col)] + resp_no_d.re, resp_no_d.im);
        }
    }
    Ok(response_matrix)
}

#[cfg(test)]
mod tests {
    use crate::{FrequencyResponse, Tf};

    use super::log_space;
    use approx::assert_abs_diff_eq;
    use num_complex::c64;

    #[test]
    fn tf_and_ss_same_freq_resp() {
        let sys_tf = (Tf::s() + 2.0) / (1.0 + Tf::s());
        let sys_ss = sys_tf.to_ss().unwrap();

        for freq in log_space(0.1, 1000.0, 10000, 10) {
            let freq_complex = c64(0.0, freq);
            let resp_tf = sys_tf.freq_response(&freq_complex);
            let resp_ss = sys_ss.freq_response(&freq_complex);
            assert_abs_diff_eq!(resp_tf.re, resp_ss.re, epsilon = 1e-3);
            assert_abs_diff_eq!(resp_tf.im, resp_ss.im, epsilon = 1e-3);
        }

        let sys_tf = sys_tf * 1.0 / Tf::s();
        let sys_ss = sys_tf.to_ss().unwrap();
        let bode_tf = sys_tf.bode(0.01, 1000.0);
        let bode_ss = sys_ss.bode(0.01, 1000.0);

        for (res_tf, res_ss) in bode_tf.iter().zip(bode_ss.iter()) {
            let mag_tf = res_tf[0];
            let phase_tf = res_tf[1];
            let mag_ss = res_ss[0];
            let phase_ss = res_ss[1];

            assert_abs_diff_eq!(mag_tf, mag_ss, epsilon = 1e-3);
            assert_abs_diff_eq!(phase_tf, phase_ss, epsilon = 1e-3);
        }

        let sys_tf = 1.0 / (0.1 + Tf::s()).powi(10);
        let sys_ss = sys_tf.to_ss().unwrap();

        let nyq_tf = sys_tf.nyquist(0.01, 1000.0);
        let nyq_ss = sys_ss.nyquist(0.01, 1000.0);

        for (res_tf, res_ss) in nyq_tf.iter().zip(nyq_ss.iter()) {
            assert_abs_diff_eq!(res_tf.re, res_ss.re, epsilon = 1e-3);
            assert_abs_diff_eq!(res_tf.im, res_tf.im, epsilon = 1e-3);
        }
    }
}
