# Control Systems TorBox

[![Build Status]][ci_action] [![Latest Version]][crates.io] [![codecov]][codecov_num]

[Latest Version]: https://img.shields.io/crates/v/control_systems_torbox.svg
[crates.io]: https://crates.io/crates/control_systems_torbox

[codecov]: https://codecov.io/gh/TorBorve/Control-Systems-TorBox/graph/badge.svg?token=VOW8UCR3FI
[codecov_num]: https://codecov.io/gh/TorBorve/Control-Systems-TorBox

[Build Status]: https://github.com/TorBorve/Control-Systems-TorBox/actions/workflows/ci.yml/badge.svg
[ci_action]: https://github.com/TorBorve/Control-Systems-TorBox/actions/workflows/ci.yml

Control Systems Torbox is a rust library for designing, analysing, simulating and implementating linear control systems.

## Usage

To use this crate add `control-systems-torbox` to your `Cargo.toml`:

```toml
[dependencies]
control-systems-torbox = "*"
```

## Examples

The crate contains a representation for transfer functions:

```rust
use control_systems_torbox::*;

let tf_lp = Tf::<f64, Continuous>::new(&[10.0], &[10.0, 1.0]); // 10/(s + 10)
let sys = 1.0 / Tf::s();
let pi_c = 1.0 + 1.0 / Tf::s();

let openloop = sys * pi_c * tf_lp;
let tracking = openloop.feedback(&Tf::new_from_scalar(1.0));
```

State-Space representations is also implemented and you can convert between transfer function and state-space

```rust
use control_systems_torbox::*;
use nalgebra::dmatrix;

let a = dmatrix![0., 1.;
                 0., 0.];
let b = dmatrix![0.;
                 1.];
let c = dmatrix![1., 0.];
let d = dmatrix![0.];

let sys_ss = Ss::<Continuous>::new(a, b, c, d).unwrap();
let sys_tf = sys_ss.to_tf().unwrap();
```

System norms such as H2-norm and H-inf-norm are also possible to calculte

```rust
use control_systems_torbox::*;

let low_pass = 1.0/(Tf::s().powi(2) + 2.0*0.2*Tf::s() + 0.1);
let low_pass = low_pass.to_ss().unwrap();
let h2_norm  = low_pass.norm_h2().unwrap();
let hinf_norm = low_pass.norm_hinf().unwrap();
println!("H2 norm is: {}, H-inf norm is: {}", h2_norm, hinf_norm);
```

It is also possible to calculte the poles and zeros of a system

```rust
use control_systems_torbox::*;

let tf = (Tf::s() -1.0)*(Tf::s() + 3.0)/((Tf::s() - 4.0)*(Tf::s() + 5.0));
let ss = tf.to_ss().unwrap();
let poles = ss.poles();
let zeros = ss.zeros().unwrap();
println!("Tf:\n{}", tf);
println!("Poles:");
for pole in poles {
    println!("{}", pole);
}
println!("\nZeros:");
for zero in zeros {
    println!("{}", zero);
}
```

Both Bodeplot and Nyquistplot is implemented and can be displayed natively:

```rust,no_run
use control_systems_torbox::*;

let sys = (1.0/Tf::s()) * (1.0 + 1.0/Tf::s());
let sys_clp =  sys.clone() / (1.0 + sys.clone());

let mut bode_plot = BodePlot::new(BodePlotOptions::default());
bode_plot.add_system(BodePlotData::new(sys_clp.bode(0.1, 10.)));
bode_plot.show(800, 600, "Bode plot demo").unwrap();

let mut nyq_plot = NyquistPlot::new(NyquistPlotOptions::default());
nyq_plot.add_system(NyquistPlotData::new(sys.nyquist(0.1, 100.)));
nyq_plot.show(600, 600, "Nyquist plot demo").unwrap();
```