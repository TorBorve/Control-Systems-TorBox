use control_systems_torbox as cst;

use cst::{plot::*, systems::Tf, utils::traits::Continuous, FrequencyResponse};

fn main() {
    let sys: Tf<f64, Continuous> = Tf::new(&[1.0], &[0.0, 1.0]);

    let pade: Tf<f64, Continuous> = Tf::new(&[1.0, -1.0], &[1.0, 1.0]);
    let pade2 = pade.clone() * Tf::new(&[1.0, -0.5], &[1.0, 0.5]);
    let sys_p1 = sys.clone() * pade;
    let sys_p2 = sys.clone() * pade2;

    let nyq_opts = NyquistPlotOptions::default().set_y_limits([-2., 2.]);
    let mut nyq_plot = NyquistPlot::new(nyq_opts);
    let nyq_data = sys.nyquist(0.01, 100.);
    nyq_plot.add_system(nyq_data.into());

    let data = NyquistPlotData::new(sys_p1.nyquist(0.01, 100.))
        .set_name("Pade 1")
        .set_color(RGBAColor::RED);
    nyq_plot.add_system(data);
    let data = NyquistPlotData::new(sys_p2.nyquist(0.01, 100.))
        .set_name("Pade 2")
        .set_color(RGBAColor::GREEN);
    nyq_plot.add_system(data);

    nyq_plot.show(800, 600, "Nyquist").unwrap();

    /////////////// BODE
    let bodeopts = BodePlotOptions::default();
    let bodeopts = bodeopts.set_x_limits([1., 100.]);
    let mut bodeplot = BodePlot::new(bodeopts);

    let sys: Tf<f64, Continuous> = Tf::new(&[0.0, 1.0], &[1.0, 1.0]);
    let mag_phase_freq_vec = sys.bode(0.1, 100.);

    bodeplot.add_system(
        BodePlotData::from(mag_phase_freq_vec.clone()).set_name("System 1"),
    );

    let data = BodePlotData::new(mag_phase_freq_vec.clone())
        .set_name("System 1 Green")
        .set_color(RGBAColor::GREEN);
    bodeplot.add_system(data);

    let sys2: Tf<f64, Continuous> = Tf::new(&[1.0], &[1.0, 1.0]);
    let mag_phase_freq_vec = sys2.bode(0.1, 100.);

    bodeplot.add_system(mag_phase_freq_vec.into());

    bodeplot.show(800, 600, "Bodeplot 1").unwrap();
}
