use control_systems_torbox as cst;
use cst::{frequency_response::*, plot::*, tf::Tf, traits::*};

fn main() {
    let sys: Tf<f64, Continuous> = Tf::new(&[1.0], &[0.0, 1.0]);

    let pade: Tf<f64, Continuous> = Tf::new(&[1.0, -1.0], &[1.0, 1.0]);
    let pade2 = pade.clone() * Tf::new(&[1.0, -0.5], &[1.0, 0.5]);
    let sys_p1 = sys.clone() * pade;
    let sys_p2 = sys.clone() * pade2;

    let nyq_opts = NyquistPlotOptions::default()
        .set_x_limits([-2., 1.])
        .set_y_limits([-2., 2.]);
    let mut nyq_plot = NyquistPlot::new(nyq_opts);
    let nyq_data = nyquist(sys, 0.01, 100.);
    nyq_plot.add_system(nyq_data.into());

    let data = NyquistPlotData::new(nyquist(sys_p1, 0.01, 100.))
        .set_name("Pade 1")
        .set_color(RGBAColor::RED);
    nyq_plot.add_system(data);
    let data = NyquistPlotData::new(nyquist(sys_p2, 0.01, 100.))
        .set_name("Pade 2")
        .set_color(RGBAColor::GREEN);
    nyq_plot.add_system(data);

    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([800., 600.]),
        ..Default::default()
    };
    let nyq_plot_egui = NyquistPlotEgui::new(nyq_plot);
    eframe::run_native(
        "Nyqust Plot",
        opts.clone(),
        Box::new(move |_cc| Ok(Box::new(nyq_plot_egui))),
    )
    .unwrap();

    let bodeopts = BodePlotOptions::default();
    let mut bodeplot = BodePlot::new(bodeopts);

    let sys: Tf<f64, Continuous> = Tf::new(&[0.0, 1.0], &[1.0, 1.0]);
    let (mag, phase, freq) = bode(sys, 0.1, 100.);
    let mag_phase_freq_vec: Vec<[f64; 3]> = mag
        .iter()
        .zip(phase.iter().zip(freq.iter()))
        .map(|(&mag, (&phase, &freq))| [mag, phase, freq])
        .collect();

    bodeplot.add_system(
        BodePlotData::from(mag_phase_freq_vec.clone()).set_name("Test"),
    );

    bodeplot.add_system(
        BodePlotData::new(mag_phase_freq_vec.clone()).set_name("INLINE"),
    );
    let data = BodePlotData::new(mag_phase_freq_vec.clone())
        .set_name("Test")
        .set_color(RGBAColor::GREEN);
    bodeplot.add_system(data);

    let mut data: BodePlotData = mag_phase_freq_vec.clone().into();
    data = data.set_color(RGBAColor::RED).set_name("TEST");
    bodeplot.add_system(data);

    bodeplot.add_system(mag_phase_freq_vec.into());
    let sys2: Tf<f64, Continuous> = Tf::new(&[1.0], &[1.0, 1.0]);
    let (mag, phase, freq) = bode(sys2, 0.1, 100.);
    let mag_phase_freq_vec: Vec<[f64; 3]> = mag
        .iter()
        .zip(phase.iter().zip(freq.iter()))
        .map(|(&mag, (&phase, &freq))| [mag, phase, freq])
        .collect();

    bodeplot.add_system(mag_phase_freq_vec.into());

    let bodeegui = BodePlotEgui::new(bodeplot);
    eframe::run_native(
        "Bode Plot",
        opts,
        Box::new(move |_cc| Ok(Box::new(bodeegui))),
    )
    .unwrap();
}
