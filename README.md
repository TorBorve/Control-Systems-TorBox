# Control Systems Torbox

Control Systems Torbox is a rust library for performing analysis, simulation and implementation of linear control systems.


## Features

- [X] Transfer functions
- [ ] State Space
- [ ] General Linear Time Invariant systems. (delays in feedback)
- [ ] Connecting and feedback of LTI systems.
- [ ] poles and zeros of LTI systems
- [X] Bodeplot
- [X] Nyquistplot
- [ ] pole zero plot
- [ ] Pade approximation
- [ ] c2d and d2c
- [ ] Simulate systems

## Todo

- [x] Display TF, RF, Polynomial
- [ ] 2.0*Tf, ...
- [ ] Improve plotting structure. Make plotting trait which nyquistplot and bodeplot implements?
- [ ] Fix pixel discretization of plotters?
- [ ] Impl egui for bode- and nyquist-plot.
- [ ] Fix x_limits and y_limits egui use auto scale if not set.

## Resources

- [plotters tutorial](https://github.com/wiseaidev/rust-data-analysis/blob/main/6-plotters-tutorial-part-1.ipynb)