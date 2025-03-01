use core::f64;

use plotters::{
    chart::{ChartBuilder, LabelAreaPosition},
    prelude::{IntoDrawingArea, *},
    series::LineSeries,
    style::WHITE,
};

fn color_order(idx: usize) -> RGBAColor {
    // Using same as MATLAB
    let colors = [
        [0., 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ];

    let idx = idx % colors.len();
    let [r, g, b] = colors[idx];
    RGBAColor((r * 255.) as u8, (g * 255.) as u8, (b * 255.) as u8, 1.0)
}

#[derive(Clone)]
pub struct BodePlotData {
    mag_phase_freq_points: Vec<[f64; 3]>,
    name: String,
    color: Option<RGBAColor>,
}

impl BodePlotData {
    pub fn new(mag_phase_freq_points: Vec<[f64; 3]>) -> Self {
        Self {
            mag_phase_freq_points,
            name: "".to_string(),
            color: None,
        }
    }

    pub fn set_name(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }
    pub fn set_color(mut self, color: RGBAColor) -> Self {
        self.color = Some(color);
        self
    }
}

impl From<Vec<[f64; 3]>> for BodePlotData {
    fn from(value: Vec<[f64; 3]>) -> Self {
        BodePlotData::new(value)
    }
}

#[derive(Default, Clone)]
pub struct BodePlotOptions {
    title: String,
    x_limits: Option<[f64; 2]>,
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

        let ([mag_min, mag_max], [phase_min, phase_max], [freq_min, freq_max]) =
            self.determine_axis_limits();

        let freq_min = freq_min.log10();
        let freq_max = freq_max.log10();

        let (mag_plot, phase_plot) =
            root.split_vertically(root.dim_in_pixel().1 / 2);

        let mut mag_chart = ChartBuilder::on(&mag_plot)
            .caption(&self.options.title, ("sans-serif", 30))
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(freq_min..freq_max, mag_min..mag_max)
            .unwrap();
        mag_chart.configure_mesh().draw().unwrap();
        for (idx, data) in self.plot_data.iter().enumerate() {
            let color = data.color.unwrap_or(color_order(idx));

            let mag_data = data
                .mag_phase_freq_points
                .iter()
                .map(|x| (x[2].log10(), x[0]));
            let series = mag_chart
                .draw_series(LineSeries::new(mag_data, color))
                .unwrap();

            if !data.name.is_empty() {
                series.label(&data.name).legend(move |(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], color)
                });
            }

            mag_chart
                .configure_series_labels()
                .background_style(&WHITE)
                .border_style(&BLACK)
                .draw()
                .unwrap();
        }

        let mut phase_chart = ChartBuilder::on(&phase_plot)
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(freq_min..freq_max, phase_min..phase_max)
            .unwrap();
        phase_chart.configure_mesh().draw().unwrap();
        for (idx, data) in self.plot_data.iter().enumerate() {
            let color = data.color.unwrap_or(color_order(idx));
            let phase_data = data
                .mag_phase_freq_points
                .iter()
                .map(|x| (x[2].log10(), x[1]));
            phase_chart
                .draw_series(LineSeries::new(phase_data, color))
                .unwrap();
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
        assert!(
            idx < 3,
            "Index out of range. Only mag, phase, freq can be indexed"
        );

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

    use super::*;
    use crate::{bode::bode, tf::Tf, traits::Continuous};
    use gtk::{cairo, prelude::*};
    use plotters::style::full_palette::PURPLE;
    use plotters_svg::SVGBackend;

    #[test]
    fn plot() {
        // let backend =plotters_cairo::CairoBackend::new();
        // let backend =
        //     plotters::backend::BitMapBackend::new("test_file.png", (2000,
        // 2000));
        let backend = SVGBackend::new("test_fig.svg", (800, 600));

        let bodeopts = BodePlotOptions::default();
        let mut bodeplot = BodePlot::new(bodeopts);

        let sys: Tf<f64, Continuous> = Tf::new(&[0.0, 1.0], &[1.0, 1.0]);
        let (mag, phase, freq) = bode(sys, 0.1, 100.);
        let mag_phase_freq_vec: Vec<[f64; 3]> = mag
            .iter()
            .zip(phase.iter().zip(freq.iter()))
            .map(|(&mag, (&phase, &freq))| [mag, phase, freq])
            .collect();

        // let mut bode_data =
        // BodePlotData::new(mag_phase_freq_vec.clone()).set_name("Test");
        // // bodebode_data.set_name("LP");

        // bodeplot.add_system(*bode_data);

        bodeplot.add_system(
            BodePlotData::from(mag_phase_freq_vec.clone()).set_name("Test"),
        );

        bodeplot.add_system(
            BodePlotData::new(mag_phase_freq_vec.clone()).set_name("INLINE"),
        );
        let data = BodePlotData::new(mag_phase_freq_vec.clone())
            .set_name("Test")
            .set_color(PURPLE.to_rgba());
        bodeplot.add_system(data);

        let mut data: BodePlotData = mag_phase_freq_vec.clone().into();
        let data = data.set_color(BLACK.into()).set_name("TEST");
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
        // let bode_data = BodePlotData {
        //     mag_phase_freq_points: mag_phase_freq_vec,
        //     name: "HP".to_string(),
        //     color: RED.to_rgba(),
        // };
        // bodeplot.add_system(bode_data);

        bodeplot.plot(backend).unwrap();

        // GTK
        let app = gtk::Application::new(Some("gtk.demo"), Default::default());

        // let bodeplot = std::rc::Rc::new(bodeplot);

        app.connect_activate(move |app| {
            let window = gtk::ApplicationWindow::new(app);
            window.set_title(Some("Bodeplot"));
            window.set_default_size(2 * 800, 2 * 600);

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
