use num_complex::{Complex64, c64};

use crate::{
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

impl<U: Time> Tf<f64, U> {
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
    pub fn bode(&self, min_freq: f64, max_freq: f64) -> Vec<[f64; 3]> {
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
    pub fn bode_freqs(
        &self,
        freqs: impl Iterator<Item = f64>,
    ) -> Vec<[f64; 3]> {
        let mut mag_phase_freq_vec = Vec::with_capacity(freqs.size_hint().0);

        for omega in freqs {
            let c = c64(0., omega);
            let sys_val = self.eval(&c);
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
    pub fn nyquist(&self, min_freq: f64, max_freq: f64) -> Vec<Complex64> {
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
    pub fn nyquist_freqs(
        &self,
        freqs: impl Iterator<Item = f64>,
    ) -> Vec<Complex64> {
        let mut pos_vals = Vec::with_capacity(freqs.size_hint().0);
        let mut neg_vals = Vec::with_capacity(freqs.size_hint().0);

        for freq in freqs {
            pos_vals.push(self.eval(&c64(0., freq)));
            neg_vals.push(self.eval(&c64(0., -freq)));
        }

        pos_vals.extend(neg_vals.iter().rev());
        pos_vals
    }
}
