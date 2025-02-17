use std::ops::Add;
use std::ops::Div;
use std::ops::Sub;
use std::vec;

use std::fmt;
use std::ops::Mul;

#[derive(Clone, PartialEq, Debug)]
pub struct Tf {
    numerator: Vec<f64>,
    denominator: Vec<f64>,
    sampling_time: f64,
}

impl Tf {
    pub fn new(
        numerator: Vec<f64>,
        denominator: Vec<f64>,
        sampling_time: f64,
    ) -> Result<Tf, String> {
        if sampling_time < 0. {
            return Err(
                "Sampling time must be zero for cotinous systems and positve for discrete systems"
                    .to_string(),
            );
        }
        if numerator.is_empty() || denominator.is_empty() {
            return Err("Numerator and Denominator cannot be empty".to_string());
        }
        if !numerator.iter().any(|&x| x != 0.) {
            return Err("Denominator cannot contain only zeros".to_string());
        }

        Ok(Tf {
            numerator,
            denominator,
            sampling_time,
        })
    }

    pub fn new_const(const_val: f64, sampling_time: f64) -> Result<Tf, String> {
        Tf::new(vec![const_val], vec![1.0], sampling_time)
    }

    pub fn s() -> Tf {
        Tf::new(vec![1.0, 0.0], vec![1.0], 0.0).unwrap()
    }

    pub fn z(sampling_time: f64) -> Result<Tf, String> {
        if sampling_time == 0.0 {
            return Err(
                "Sampling time must be positive for discrete transfer function".to_string(),
            );
        }
        Tf::new(vec![1.0, 0.0], vec![1.0], sampling_time)
    }

    pub fn inv(self) -> Self {
        Tf {
            numerator: self.denominator,
            denominator: self.numerator,
            sampling_time: self.sampling_time,
        }
    }
}

impl Mul for Tf {
    type Output = Result<Self, String>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.sampling_time != rhs.sampling_time {
            return Err("Sampling times must agree".to_string());
        }
        let new_numerator = mul_polynomials(&self.numerator, &rhs.numerator);
        let new_denom = mul_polynomials(&self.denominator, &rhs.denominator);

        Tf::new(new_numerator, new_denom, self.sampling_time)
    }
}

impl Div for Tf {
    type Output = Result<Self, String>;

    fn div(self, rhs: Self) -> Self::Output {
        let rhs_inv = rhs.inv();
        self * rhs_inv
    }
}

impl Add for Tf {
    type Output = Result<Self, String>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.sampling_time != rhs.sampling_time {
            return Err("Sampling times must agree".to_string());
        }
        let common_den = mul_polynomials(&self.denominator, &rhs.denominator);
        let lhs_num = mul_polynomials(&self.numerator, &rhs.denominator);
        let rhs_num = mul_polynomials(&rhs.numerator, &self.denominator);
        let sum_num = add_polynomials(&lhs_num, &rhs_num);
        Tf::new(sum_num, common_den, self.sampling_time)
    }
}

impl Sub for Tf {
    type Output = Result<Self, String>;

    fn sub(self, rhs: Self) -> Self::Output {
        let neg_gain = Tf::new_const(-1.0, self.sampling_time)?;
        let neg_rhs = (neg_gain * rhs)?;
        self + neg_rhs
    }
}

fn mul_polynomials(lhs: &Vec<f64>, rhs: &Vec<f64>) -> Vec<f64> {
    let mut result = vec![0.; lhs.len() + rhs.len() - 1];

    for (idx_l, val_l) in lhs.iter().enumerate() {
        for (idx_r, val_r) in rhs.iter().enumerate() {
            let new_val = val_l * val_r;
            result[idx_l + idx_r] += new_val;
        }
    }

    result
}

fn add_polynomials(lhs: &Vec<f64>, rhs: &Vec<f64>) -> Vec<f64> {
    let mut result = vec![0.; lhs.len().max(rhs.len())];
    for (idx, val_l) in lhs.iter().enumerate() {
        result[idx] += val_l;
    }
    for (idx, val_r) in rhs.iter().enumerate() {
        result[idx] += val_r;
    }
    result
}

impl fmt::Display for Tf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let domain_char = if self.sampling_time == 0. { 's' } else { 'z' };
        let num_str = polynomial_to_string(&self.numerator, domain_char);
        let den_str = polynomial_to_string(&self.denominator, domain_char);

        let h_bar = "\u{2500}".repeat(num_str.len().max(den_str.len()));
        write!(f, "{}\n{}\n{}", num_str, h_bar, den_str)
    }
}

fn polynomial_to_string(poly: &Vec<f64>, char_variable: char) -> String {
    let mut poly_str = "".to_string();
    for (idx, val) in poly.iter().enumerate() {
        let order = poly.len() - 1 - idx;
        if !(*val != 0. || (order == 0 && poly_str.is_empty())) {
            continue;
        }
        let variable_str = if order > 1 {
            format!("{}^{}", char_variable, order)
        } else if order == 1 {
            char_variable.to_string()
        } else {
            "".to_string()
        };
        let sign_str = if idx == 0 {
            let mut s_str = "".to_string();
            if val.is_sign_negative() {
                s_str.push('-');
            }
            s_str
        } else {
            let s_str = if val.is_sign_negative() { " - " } else { " + " };
            s_str.to_string()
        };
        let val_str = if val.abs() == 1.0 && order != 0 {
            "".to_string()
        } else {
            val.abs().to_string()
        };
        poly_str.push_str(&format!("{}{}{}", sign_str, val_str, variable_str));
    }

    poly_str
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_tf() {
        {
            let s = Tf::new(vec![1.0, 0.0], vec![1.0], 0.0);
            assert!(s.is_ok());
            let s = s.unwrap();
            let s_new = Tf::s();
            assert_eq!(s, s_new);
        }
        {
            let t_s = 0.1;
            let z = Tf::new(vec![1.0, 0.0], vec![1.], t_s);
            assert!(z.is_ok());
            let z = z.unwrap();
            let z_new = Tf::z(0.0);
            assert!(z_new.is_err());
            let z_new = Tf::z(t_s);
            assert!(z_new.is_ok());
            let z_new = z_new.unwrap();
            assert_eq!(z, z_new);
        }
    }

    #[test]
    fn tf_operators() {
        {
            let s = Tf::s();
            let s_inv = s.clone().inv();
            let s_on_s = s * s_inv;
            let expected_res = Tf::new(vec![1., 0.], vec![1., 0.], 0.);
            assert_eq!(s_on_s, expected_res);
        }
        {
            let one = Tf::new_const(1.0, 0.0).unwrap();
            let two = Tf::new_const(2.0, 0.0).unwrap();
            assert_eq!((one.clone() + one.clone()).unwrap(), two);
            assert_eq!((two.clone() - one.clone()).unwrap(), one);
        }
        {}
    }
}
