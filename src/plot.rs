use core::f64;

use egui::Color32;
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotPoints};
use num_complex::Complex64;
use plotters::{
    chart::{ChartBuilder, LabelAreaPosition},
    prelude::{IntoDrawingArea, *},
    series::LineSeries,
    style::WHITE,
};

// See: https://github.com/wiseaidev/rust-data-analysis/blob/main/6-plotters-tutorial-part-1.ipynb

#[allow(clippy::approx_constant)]
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

impl BodePlotOptions {
    pub fn new(title: &str) -> Self {
        Self {
            title: title.to_string(),
            x_limits: None,
        }
    }

    pub fn set_title(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }

    pub fn set_x_limits(mut self, x_limits: [f64; 2]) -> Self {
        self.x_limits = Some(x_limits);
        self
    }
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

        // let freq_min = freq_min.log10();
        // let freq_max = freq_max.log10();

        let (mag_plot, phase_plot) =
            root.split_vertically(root.dim_in_pixel().1 / 2);

        let mut mag_chart = ChartBuilder::on(&mag_plot)
            .caption(&self.options.title, ("sans-serif", 30))
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(
                (freq_min..freq_max).log_scale(),
                mag_min..mag_max,
            )
            .unwrap();
        mag_chart.configure_mesh().draw().unwrap();
        for (idx, data) in self.plot_data.iter().enumerate() {
            let color = data.color.unwrap_or(color_order(idx));

            let mag_data =
                data.mag_phase_freq_points.iter().map(|x| (x[2], x[0]));
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
                .background_style(WHITE)
                .border_style(BLACK)
                .draw()
                .unwrap();
        }

        let mut phase_chart = ChartBuilder::on(&phase_plot)
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(
                (freq_min..freq_max).log_scale(),
                phase_min..phase_max,
            )
            .unwrap();
        phase_chart.configure_mesh().draw().unwrap();
        for (idx, data) in self.plot_data.iter().enumerate() {
            let color = data.color.unwrap_or(color_order(idx));
            let phase_data =
                data.mag_phase_freq_points.iter().map(|x| (x[2], x[1]));
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

#[derive(Clone, Debug)]
pub struct NyquistPlotData {
    data_points: Vec<Complex64>,
    name: String,
    color: Option<RGBAColor>,
}

impl NyquistPlotData {
    pub fn new(data_points: Vec<Complex64>) -> Self {
        Self {
            data_points,
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

impl From<Vec<Complex64>> for NyquistPlotData {
    fn from(value: Vec<Complex64>) -> Self {
        NyquistPlotData::new(value)
    }
}

#[derive(Default, Clone, Debug)]
pub struct NyquistPlotOptions {
    title: String,
    x_limits: Option<[f64; 2]>,
    y_limits: Option<[f64; 2]>,
}

impl NyquistPlotOptions {
    pub fn set_title(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }

    pub fn set_x_limits(mut self, x_limits: [f64; 2]) -> Self {
        self.x_limits = Some(x_limits);
        self
    }

    pub fn set_y_limits(mut self, y_limits: [f64; 2]) -> Self {
        self.y_limits = Some(y_limits);
        self
    }
}

#[derive(Debug, Default, Clone)]
pub struct NyquistPlot {
    options: NyquistPlotOptions,
    plot_data: Vec<NyquistPlotData>,
}

impl NyquistPlot {
    pub fn new(options: NyquistPlotOptions) -> Self {
        Self {
            options,
            plot_data: vec![],
        }
    }

    pub fn add_system(&mut self, data: NyquistPlotData) -> &mut Self {
        self.plot_data.push(data);
        self
    }

    pub fn plot<B: DrawingBackend>(
        &self,
        backend: B,
    ) -> Result<(), Box<dyn std::error::Error>>
    where
        B::ErrorType: std::error::Error + 'static,
    {
        let ([x_min, x_max], [y_min, y_max]) =
            NyquistPlot::determine_axis_limits(
                self.plot_data
                    .iter()
                    .flat_map(|data| data.data_points.iter()),
                self.options.x_limits,
                self.options.y_limits,
            );

        let root = backend.into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption(&self.options.title, ("sans-serif", 30))
            .margin(20)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

        chart.configure_mesh().draw()?;
        for (idx, data) in self.plot_data.iter().enumerate() {
            let color = data.color.unwrap_or(color_order(idx));

            let plot_data = data.data_points.iter().map(|c| (c.re, c.im));
            let series =
                chart.draw_series(LineSeries::new(plot_data, color))?;

            if !data.name.is_empty() {
                series.label(&data.name).legend(move |(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], color)
                });
            }
        }

        chart
            .configure_series_labels()
            .background_style(WHITE)
            .border_style(BLACK)
            .draw()?;

        // TODO: Draw marker for (-1, 0)
        Ok(())
    }

    fn determine_axis_limits<'a, I>(
        data: I,
        x_limits: Option<[f64; 2]>,
        y_limits: Option<[f64; 2]>,
    ) -> ([f64; 2], [f64; 2])
    where
        I: Iterator<Item = &'a Complex64>,
    {
        let mut x_min_auto = f64::MAX;
        let mut x_max_auto = f64::MIN;
        let mut y_min_auto = f64::MAX;
        let mut y_max_auto = f64::MIN;

        for c in data {
            x_min_auto = x_min_auto.min(c.re);
            x_max_auto = x_max_auto.max(c.re);
            y_min_auto = y_min_auto.min(c.im);
            y_max_auto = y_max_auto.max(c.im);
        }

        let x_limits = x_limits.unwrap_or([x_min_auto, x_max_auto]);
        let y_limits = y_limits.unwrap_or([y_min_auto, y_max_auto]);
        (x_limits, y_limits)
    }
}

pub struct NyquistPlotEgui {
    nyq: NyquistPlot,
    init: bool,
}

impl NyquistPlotEgui {
    pub fn new(nyq: NyquistPlot) -> Self {
        Self { nyq, init: false }
    }
}

impl From<NyquistPlot> for NyquistPlotEgui {
    fn from(value: NyquistPlot) -> Self {
        Self::new(value)
    }
}

impl eframe::App for NyquistPlotEgui {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let plot = Plot::new("Nyquist Plot").legend(Legend::default());
            plot.show(ui, |plot_ui| {
                if !self.init {
                    let mut bounds = plot_ui.plot_bounds();
                    if let Some([x_min, x_max]) = self.nyq.options.x_limits {
                        bounds.set_x(&PlotBounds::from_min_max(
                            [x_min, x_min],
                            [x_max, x_max],
                        ));
                    }
                    if let Some([y_min, y_max]) = self.nyq.options.y_limits {
                        bounds.set_y(&PlotBounds::from_min_max(
                            [y_min, y_min],
                            [y_max, y_max],
                        ));
                    }
                    plot_ui.set_plot_bounds(bounds);
                }
                for data in &self.nyq.plot_data {
                    let plot_points = PlotPoints::from_iter(
                        data.data_points.iter().map(|c| [c.re, c.im]),
                    );
                    let mut line = Line::new(plot_points);
                    if !data.name.is_empty() {
                        line = line.name(data.name.clone());
                    }
                    if let Some(color) = data.color {
                        let color = Color32::from_rgba_premultiplied(
                            color.0,
                            color.1,
                            color.2,
                            (color.3 * 255.) as u8,
                        );
                        line = line.color(color);
                    }
                    plot_ui.line(line);
                }
            });
        });
        self.init = true;
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::{
        frequency_response::{bode, nyquist},
        tf::Tf,
        traits::Continuous,
    };
    use gtk::{cairo, prelude::*};
    use plotters::style::full_palette::PURPLE;
    use plotters_svg::SVGBackend;

    #[test]
    fn bodeplot() {
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
        data = data.set_color(BLACK.into()).set_name("TEST");
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

        // app.run();
    }

    #[test]
    fn nyquistplot() {
        let sys: Tf<f64, Continuous> = Tf::new(&[1.0], &[0.0, 1.0]);

        let pade: Tf<f64, Continuous> = Tf::new(&[1.0, -1.0], &[1.0, 1.0]);
        let pade2 = pade.clone() * Tf::new(&[1.0, -0.5], &[1.0, 0.5]);
        let sys_p1 = sys.clone() * pade;
        let sys_p2 = sys.clone() * pade2;

        let backend = SVGBackend::new("nyquist.svg", (800, 600));

        let mut nyq_opts = NyquistPlotOptions::default();
        nyq_opts = nyq_opts.set_x_limits([-5., 1.]).set_y_limits([-10., 10.]);
        let mut nyq_plot = NyquistPlot::new(nyq_opts);
        let nyq_data = nyquist(sys, 0.01, 100.);
        nyq_plot.add_system(nyq_data.into());

        let data = NyquistPlotData::new(nyquist(sys_p1, 0.01, 100.))
            .set_name("Pade 1");
        nyq_plot.add_system(data);
        let data = NyquistPlotData::new(nyquist(sys_p2, 0.01, 100.))
            .set_name("Pade 2");
        nyq_plot.add_system(data);

        nyq_plot.plot(backend).unwrap();

        // GTK
        // let app = gtk::Application::new(Some("gtk.demo"),
        // Default::default());

        // // let bodeplot = std::rc::Rc::new(bodeplot);

        // app.connect_activate(move |app| {
        //     let window = gtk::ApplicationWindow::new(app);
        //     window.set_title(Some("Bodeplot"));
        //     window.set_default_size(800, 600);

        //     let draw_area = gtk::DrawingArea::new();
        //     let nyq_plot = nyq_plot.clone();
        //     draw_area.set_draw_func(move |_, cr: &cairo::Context, _, _| {
        //         let backend =
        //             plotters_cairo::CairoBackend::new(cr, (800,
        // 600)).unwrap();         nyq_plot.plot(backend).unwrap();
        //     });

        //     window.set_child(Some(&draw_area));
        //     window.show();
        // });

        // app.run();

        // let opts = eframe::NativeOptions {
        //     viewport: egui::ViewportBuilder::default().with_inner_size([350.,
        // 200.]),     ..Default::default()
        // };
        // let nyq_plot_clone = nyq_plot.clone();
        // eframe::run_native("Nyqust Plot", opts, Box::new(move |_cc|
        // Ok(Box::new(nyq_plot_clone)) )).unwrap();

        let backend = SVGBackend::new("example.svg", (800, 600));

        let root = backend.into_drawing_area();
        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .caption("10 * x", ("sans-serif", 30))
            .margin(20)
            .build_cartesian_2d(-10f64..10f64, -10f64..10f64)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        let data: Vec<(f64, f64)> = (-1000..=1000)
            .map(|x| x as f64 / 1000.0)
            .map(|x| (x, 10.0 * x))
            .collect();

        chart.draw_series(LineSeries::new(data, &RED)).unwrap();

        root.present().unwrap();
    }
}
