use core::f64;

use gtk::{cairo, prelude::*};
use plotters::{
    backend,
    chart::{ChartBuilder, LabelAreaPosition},
    prelude::{IntoDrawingArea, *},
    series::LineSeries,
    style::{Color, BLUE, WHITE},
};
use plotters_backend::BackendColor;
use plotters_cairo::CairoBackend;

use crate::{tf::Tf, traits::Continuous};

#[derive(Clone)]
pub struct BodePlotData {
    mag_phase_freq_points: Vec<[f64; 3]>,
    name: String,
    color: RGBAColor,
}

#[derive(Default, Clone)]
pub struct BodePlotOptions {
    title: String,
    x_limits: Option<[f64; 2]>,
    // y_lims: Option<[f64; 2]>,
}

#[derive(Default, Clone)]
pub struct BodePlot {
    options: BodePlotOptions,
    plot_data: Vec<BodePlotData>,
}

impl BodePlot {
    pub fn new(options: BodePlotOptions) -> Self {
        Self {
            options,
            plot_data: vec![],
        }
    }

    pub fn add_system(&mut self, data: BodePlotData) -> &mut Self {
        self.plot_data.push(data);
        self
    }

    pub fn plot<B: DrawingBackend>(
        &self,
        backend: B,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let root = backend.into_drawing_area();
        root.fill(&WHITE).unwrap();

        // let x_min = self.options.x_lims.map_or(f64::MAX, |lims| lims[0]);
        let ([mag_min, mag_max], [phase_min, phase_max], [freq_min, freq_max]) = self.determine_axis_limits();

        let freq_min = freq_min.log10();
        let freq_max = freq_max.log10();

        // println!("y-min: {}, y-max: {}", y_min, y_max);

        let mut chart = ChartBuilder::on(&root)
            .caption(&self.options.title, ("sans-serif", 30))
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(freq_min..freq_max, mag_min..mag_max)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        for data in &self.plot_data {
            let mag_data = data.mag_phase_freq_points.iter().map(|x| (x[2].log10(), x[0]));
            chart.draw_series(LineSeries::new(
                mag_data,
                data.color
            )).unwrap();

        }
        Ok(())
    }

    fn determine_axis_limits(&self) -> ([f64; 2], [f64; 2], [f64; 2]) {

        let mag_limits_auto = self.determine_axis_limit_idx(0);
        let phase_limits_auto = self.determine_axis_limit_idx(1);
        let freq_limits_auto = self.determine_axis_limit_idx(2);

        let freq_limits = self.options.x_limits.unwrap_or(freq_limits_auto);

        (mag_limits_auto, phase_limits_auto, freq_limits)
    }

    fn determine_axis_limit_idx(&self, idx: usize) -> [f64; 2] {
        assert!(idx < 3, "Index out of range. Only mag, phase, freq can be indexed");

        let mut min_val = f64::MAX;
        let mut max_val = f64::MIN;

        for data in &self.plot_data {
            for p in &data.mag_phase_freq_points {
                min_val = min_val.min(p[idx]);
                max_val = max_val.max(p[idx]);
            }
        }

        [min_val, max_val]
    }



}

#[cfg(test)]
mod tests {

    use crate::bode::bode;

    use super::*;

    #[test]
    fn plot() {
        // let backend =plotters_cairo::CairoBackend::new();
        let backend =
            plotters::backend::BitMapBackend::new("test_file.png", (800, 600));

        let bodeopts = BodePlotOptions::default();
        let mut bodeplot = BodePlot::new(bodeopts);

        let sys: Tf<f64, Continuous> = Tf::new(&[0.0, 1.0], &[1.0, 1.0]);
        let (mag, phase, freq) = bode(sys, 0.1, 100.);
        let mag_phase_freq_vec: Vec<[f64; 3]> = mag.iter().zip(phase.iter().zip(freq.iter())).map(|(&mag, (&phase, &freq))| [mag, phase, freq]).collect();
        let bode_data = BodePlotData {
            mag_phase_freq_points: mag_phase_freq_vec,
            name: "System".to_string(),
            color: BLUE.to_rgba(),
        };
        bodeplot.add_system(bode_data);
        let sys2: Tf<f64, Continuous> = Tf::new(&[1.0], &[1.0, 1.0]);
        let (mag, phase, freq) = bode(sys2, 0.1, 100.);
        let mag_phase_freq_vec: Vec<[f64; 3]> = mag.iter().zip(phase.iter().zip(freq.iter())).map(|(&mag, (&phase, &freq))| [mag, phase, freq]).collect();
        let bode_data = BodePlotData {
            mag_phase_freq_points: mag_phase_freq_vec,
            name: "System".to_string(),
            color: RED.to_rgba(),
        };
        bodeplot.add_system(bode_data);


        bodeplot.plot(backend).unwrap();

        // GTK
        let app = gtk::Application::new(Some("gtk.demo"), Default::default());

        // let bodeplot = std::rc::Rc::new(bodeplot);

        app.connect_activate(move |app| {
            let window = gtk::ApplicationWindow::new(app);
            window.set_title(Some("Bodeplot"));
            window.set_default_size(800, 600);

            let draw_area = gtk::DrawingArea::new();
            let bodeplot = bodeplot.clone();
            draw_area.set_draw_func(move |_, cr: &cairo::Context, _, _| {
                let backend =
                    plotters_cairo::CairoBackend::new(cr, (800, 600)).unwrap();
                bodeplot.plot(backend).unwrap();
            });

            window.set_child(Some(&draw_area));
            window.show();
        });

        app.run();
    }
}
