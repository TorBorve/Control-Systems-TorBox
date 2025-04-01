/*!
# Control Systems Torbox

Control Systems Torbox is a Rust library for designing, analyzing, simulating, and implementing linear control systems. This crate provides a comprehensive set of tools for working with transfer functions, state-space representations, frequency response analysis, and plotting.

## Features

- **Transfer Functions**: Create and manipulate transfer functions.
- **State-Space Representations**: Convert between transfer functions and state-space representations.
- **Frequency Response Analysis**: Compute Bode and Nyquist plots.
- **Properties of System**: H2-norm, H-inf-norm, poles, zeros
- **Plotting**: Display Bode and Nyquist plots natively using `egui`.


## Examples

The crate contains a representation for transfer functions:

```rust
use control_systems_torbox::*;

let tf_lp = Tf::<f64, Continuous>::new(&[10.0], &[10.0, 1.0]); // 10/(s + 10)
let sys = 1.0 / Tf::s();
let pi_c = 1.0 + 1.0 / Tf::s();

let openloop = sys * pi_c * tf_lp;
let tracking = openloop.clone() / (1.0 + openloop);
```

State-Space representations is also implemented and you can convert between transfer function and state-space

```rust
use control_systems_torbox::*;
use nalgebra::DMatrix;

let a = DMatrix::from_row_slice(2, 2, &[0., 0., 1., 0.]);
let b = DMatrix::from_row_slice(2, 1, &[1., 0.]);
let c = DMatrix::from_row_slice(1, 2, &[0., 1.]);
let d = DMatrix::from_row_slice(1, 1, &[0.]);

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

!*/

pub mod analysis;
pub mod plot;
pub mod slicot_wrapper;
pub mod systems;
pub mod transformations;
pub mod utils;

pub use plot::{
    BodePlot, BodePlotData, BodePlotEgui, BodePlotOptions, NyquistPlot,
    NyquistPlotData, NyquistPlotEgui, NyquistPlotOptions, Plot, RGBAColor,
};
pub use systems::{Ss, Tf};
pub use transformations::SsRealization;
pub use utils::traits::{Continuous, Discrete};
