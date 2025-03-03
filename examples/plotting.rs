use control_systems_torbox as cst;
use cst::{frequency_response::*, plot::*, tf::Tf, traits::*};
use plotters::style::*;

fn main() {
    let sys: Tf<f64, Continuous> = Tf::new(&[1.0], &[0.0, 1.0]);

    let pade: Tf<f64, Continuous> = Tf::new(&[1.0, -1.0], &[1.0, 1.0]);
    let pade2 = pade.clone() * Tf::new(&[1.0, -0.5], &[1.0, 0.5]);
    let sys_p1 = sys.clone() * pade;
    let sys_p2 = sys.clone() * pade2;

    // let backend = SVGBackend::new("nyquist.svg", (800, 600));

    let nyq_opts = NyquistPlotOptions::default()
        .set_x_limits([-2., 1.])
        .set_y_limits([-2., 2.]);
    // nyq_opts = nyq_opts.set_x_limits([-5., 1.]).set_y_limits([-10., 10.]);
    let mut nyq_plot = NyquistPlot::new(nyq_opts);
    let nyq_data = nyquist(sys, 0.01, 100.);
    nyq_plot.add_system(nyq_data.into());

    let data = NyquistPlotData::new(nyquist(sys_p1, 0.01, 100.))
        .set_name("Pade 1")
        .set_color(RED.into());
    nyq_plot.add_system(data);
    let data = NyquistPlotData::new(nyquist(sys_p2, 0.01, 100.))
        .set_name("Pade 2")
        .set_color(GREEN.into());
    nyq_plot.add_system(data);

    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([800., 600.]),
        ..Default::default()
    };
    let nyq_plot_egui = NyquistPlotEgui::new(nyq_plot);
    eframe::run_native(
        "Nyqust Plot",
        opts,
        Box::new(move |_cc| Ok(Box::new(nyq_plot_egui))),
    )
    .unwrap();
}
