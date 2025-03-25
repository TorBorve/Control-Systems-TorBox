# Control Systems Torbox

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

```rust,no_run
use control_systems_torbox::*;

let tf_lp = Tf::<f64, Continuous>::new(&[10.0], &[10.0, 1.0]); // 10/(s + 10)
let sys = 1.0 / Tf::s();
let pi_c = 1.0 + 1.0 / Tf::s();

let openloop = sys * pi_c * tf_lp;
let tracking = openloop.clone() / (1.0 + openloop);
```

State-Space representations is also implemented and you can convert between transfer function and state-space

```rust,no_run
use control_systems_torbox::*;
use nalgebra::DMatrix;

let a = DMatrix::from_row_slice(2, 2, &[0., 0., 1., 0.]);
let b = DMatrix::from_row_slice(2, 1, &[1., 0.]);
let c = DMatrix::from_row_slice(1, 2, &[0., 1.]);
let d = DMatrix::from_row_slice(1, 1, &[0.]);

let sys_ss = Ss::<Continuous>::new(a, b, c, d).unwrap();
let sys_tf = ss2tf(&sys_ss).unwrap();
```

Both Bodeplot and Nyquistplot is implemented and can be displayed natively:

```rust,no_run
use control_systems_torbox::*;

let sys = (1.0/Tf::s()) * (1.0 + 1.0/Tf::s());
let sys_clp =  sys.clone() / (1.0 + sys.clone());

let mut bode_plot = BodePlot::new(BodePlotOptions::default());
bode_plot.add_system(BodePlotData::new(bode(sys_clp, 0.1, 10.)));
bode_plot.show(800, 600, "Bode plot demo").unwrap();

let mut nyq_plot = NyquistPlot::new(NyquistPlotOptions::default());
nyq_plot.add_system(NyquistPlotData::new(nyquist(sys, 0.1, 100.)));
nyq_plot.show(600, 600, "Nyquist plot demo").unwrap();
```
