use core::f64;

use egui::Color32;
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotPoints};
use num_complex::Complex64;

// See: https://github.com/wiseaidev/rust-data-analysis/blob/main/6-plotters-tutorial-part-1.ipynb

#[derive(Clone, Copy, Debug)]
pub struct RGBAColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}


impl RGBAColor {
    // Common colors
    pub const WHITE: RGBAColor = RGBAColor { r: 255, g: 255, b: 255, a: 255 };
    pub const BLACK: RGBAColor = RGBAColor { r: 0, g: 0, b: 0, a: 255 };
    pub const BLUE: RGBAColor = RGBAColor { r: 0, g: 0, b: 255, a: 255 };
    pub const RED: RGBAColor = RGBAColor { r: 255, g: 0, b: 0, a: 255 };
    pub const GREEN: RGBAColor = RGBAColor { r: 0, g: 255, b: 0, a: 255 };
    pub const YELLOW: RGBAColor = RGBAColor { r: 255, g: 255, b: 0, a: 255 };
    pub const CYAN: RGBAColor = RGBAColor { r: 0, g: 255, b: 255, a: 255 };
    pub const MAGENTA: RGBAColor = RGBAColor { r: 255, g: 0, b: 255, a: 255 };
    pub const ORANGE: RGBAColor = RGBAColor { r: 255, g: 165, b: 0, a: 255 };
    pub const PURPLE: RGBAColor = RGBAColor { r: 128, g: 0, b: 128, a: 255 };
    pub const BROWN: RGBAColor = RGBAColor { r: 165, g: 42, b: 42, a: 255 };
    pub const PINK: RGBAColor = RGBAColor { r: 255, g: 192, b: 203, a: 255 };
    pub const GRAY: RGBAColor = RGBAColor { r: 128, g: 128, b: 128, a: 255 };
    pub const LIGHT_GRAY: RGBAColor = RGBAColor { r: 211, g: 211, b: 211, a: 255 };
    pub const DARK_GRAY: RGBAColor = RGBAColor { r: 169, g: 169, b: 169, a: 255 };


    pub fn new(red: u8, green: u8, blue: u8, alpha: u8) -> Self {
        Self {
            r: red,
            g: green,
            b: blue,
            a: alpha,
        }
    }

    fn to_egui(&self) -> Color32 {
        Color32::from_rgba_premultiplied(self.r, self.g, self.b, self.a)
    }
}

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
    RGBAColor::new((r * 255.) as u8, (g * 255.) as u8, (b * 255.) as u8, 255)
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
}

#[derive(Clone)]
pub struct BodePlotEgui {
    bode: BodePlot,
    init: bool,
}

impl BodePlotEgui {
    pub fn new(bode: BodePlot) -> Self {
        Self { bode, init: false }
    }
}

impl From<BodePlot> for BodePlotEgui {
    fn from(value: BodePlot) -> Self {
        Self::new(value)
    }
}

impl eframe::App for BodePlotEgui {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let plot_height = ui.available_height() / 2.0;
            // let mut x_bounds = None;
            ui.vertical(|ui| {
                ui.add_sized(
                    [ui.available_width(), plot_height],
                    |ui: &mut egui::Ui| -> egui::Response {
                        let mag_plot = Plot::new("Magnitude Plot")
                            .legend(Legend::default());
                        let resp = mag_plot.show(ui, |plot_ui| {
                            if !self.init {
                                let mut bounds = plot_ui.plot_bounds();
                                if let Some([x_min, x_max]) =
                                    self.bode.options.x_limits
                                {
                                    bounds.set_x(&PlotBounds::from_min_max(
                                        [x_min, x_min],
                                        [x_max, x_max],
                                    ));
                                }
                                plot_ui.set_plot_bounds(bounds);
                            }

                            for data in &self.bode.plot_data {
                                let plot_points = PlotPoints::from_iter(
                                    data.mag_phase_freq_points
                                        .iter()
                                        .map(|p| [p[2].log10(), p[0]]),
                                );
                                let mut line = Line::new(plot_points);
                                if !data.name.is_empty() {
                                    line = line.name(data.name.clone());
                                }
                                if let Some(color) = data.color {
                                    line = line.color(color.to_egui());
                                }
                                plot_ui.line(line);
                            }
                        });
                        resp.response
                    },
                );

                ui.add_sized(
                    [ui.available_width(), plot_height],
                    |ui: &mut egui::Ui| -> egui::Response {
                        let mag_plot =
                            Plot::new("Phase Plot").legend(Legend::default());
                        mag_plot
                            .show(ui, |plot_ui| {
                                if !self.init {
                                    let mut bounds = plot_ui.plot_bounds();
                                    if let Some([x_min, x_max]) =
                                        self.bode.options.x_limits
                                    {
                                        bounds.set_x(
                                            &PlotBounds::from_min_max(
                                                [x_min, x_min],
                                                [x_max, x_max],
                                            ),
                                        );
                                    }
                                    plot_ui.set_plot_bounds(bounds);
                                }
                                for data in &self.bode.plot_data {
                                    let plot_points = PlotPoints::from_iter(
                                        data.mag_phase_freq_points
                                            .iter()
                                            .map(|p| [p[2].log10(), p[1]]),
                                    );
                                    let mut line = Line::new(plot_points);
                                    if !data.name.is_empty() {
                                        line = line.name(data.name.clone());
                                    }
                                    if let Some(color) = data.color {
                                        line = line.color(color.to_egui());
                                    }
                                    plot_ui.line(line);
                                }
                            })
                            .response
                    },
                );
            });
        });
        self.init = true;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// Nyquist
//////////////////////////////////////////////////////////////////////////////////////////////
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
                        line = line.color(color.to_egui());
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

    #[test]
    fn bodeplot() {
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
            .set_color(RGBAColor::PURPLE);
        bodeplot.add_system(data);

        let mut data: BodePlotData = mag_phase_freq_vec.clone().into();
        data = data.set_color(RGBAColor::BLACK).set_name("TEST");
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
    }

    #[test]
    fn nyquistplot() {
        let sys: Tf<f64, Continuous> = Tf::new(&[1.0], &[0.0, 1.0]);

        let pade: Tf<f64, Continuous> = Tf::new(&[1.0, -1.0], &[1.0, 1.0]);
        let pade2 = pade.clone() * Tf::new(&[1.0, -0.5], &[1.0, 0.5]);
        let sys_p1 = sys.clone() * pade;
        let sys_p2 = sys.clone() * pade2;

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
    }
}
