use core::f64;

use egui::{Align, Color32, Layout, Ui};
use egui_plot::{Legend, Line, PlotBounds, PlotPoints, PlotUi};
use num_complex::Complex64;

#[derive(Clone, Copy, Debug)]
pub struct RGBAColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl RGBAColor {
    // Common colors
    pub const WHITE: RGBAColor = RGBAColor {
        r: 255,
        g: 255,
        b: 255,
        a: 255,
    };
    pub const BLACK: RGBAColor = RGBAColor {
        r: 0,
        g: 0,
        b: 0,
        a: 255,
    };
    pub const BLUE: RGBAColor = RGBAColor {
        r: 0,
        g: 0,
        b: 255,
        a: 255,
    };
    pub const RED: RGBAColor = RGBAColor {
        r: 255,
        g: 0,
        b: 0,
        a: 255,
    };
    pub const GREEN: RGBAColor = RGBAColor {
        r: 0,
        g: 255,
        b: 0,
        a: 255,
    };
    pub const YELLOW: RGBAColor = RGBAColor {
        r: 255,
        g: 255,
        b: 0,
        a: 255,
    };
    pub const CYAN: RGBAColor = RGBAColor {
        r: 0,
        g: 255,
        b: 255,
        a: 255,
    };
    pub const MAGENTA: RGBAColor = RGBAColor {
        r: 255,
        g: 0,
        b: 255,
        a: 255,
    };
    pub const ORANGE: RGBAColor = RGBAColor {
        r: 255,
        g: 165,
        b: 0,
        a: 255,
    };
    pub const PURPLE: RGBAColor = RGBAColor {
        r: 128,
        g: 0,
        b: 128,
        a: 255,
    };
    pub const BROWN: RGBAColor = RGBAColor {
        r: 165,
        g: 42,
        b: 42,
        a: 255,
    };
    pub const PINK: RGBAColor = RGBAColor {
        r: 255,
        g: 192,
        b: 203,
        a: 255,
    };
    pub const GRAY: RGBAColor = RGBAColor {
        r: 128,
        g: 128,
        b: 128,
        a: 255,
    };
    pub const LIGHT_GRAY: RGBAColor = RGBAColor {
        r: 211,
        g: 211,
        b: 211,
        a: 255,
    };
    pub const DARK_GRAY: RGBAColor = RGBAColor {
        r: 169,
        g: 169,
        b: 169,
        a: 255,
    };

    pub fn new(red: u8, green: u8, blue: u8, alpha: u8) -> Self {
        Self {
            r: red,
            g: green,
            b: blue,
            a: alpha,
        }
    }

    fn to_egui(self) -> Color32 {
        Color32::from_rgba_premultiplied(self.r, self.g, self.b, self.a)
    }
}

#[allow(clippy::approx_constant, unused)]
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

pub trait Plot {
    fn show(
        self,
        width: usize,
        height: usize,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error + 'static>>;
}

impl<T: eframe::App> Plot for T {
    fn show(
        self,
        width: usize,
        height: usize,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let opts = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default()
                .with_inner_size([width as f32, height as f32]),
            ..Default::default()
        };
        eframe::run_native(
            name,
            opts,
            Box::new(move |_cc| Ok(Box::new(self))),
        )?;
        Ok(())
    }
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
    shared_x_limits: [f64; 2],
    last_mag_x_limits: [f64; 2],
    last_phase_x_limits: [f64; 2],
}

impl BodePlotEgui {
    pub fn new(bode: BodePlot) -> Self {
        let bounds = BodePlotEgui::get_mag_bounds(
            bode.options.x_limits,
            None,
            bode.plot_data
                .iter()
                .flat_map(|d| d.mag_phase_freq_points.iter()),
        );
        let x_limits = bounds_to_x_limits(&bounds);
        Self {
            bode,
            init: false,
            shared_x_limits: x_limits,
            last_mag_x_limits: x_limits,
            last_phase_x_limits: x_limits,
        }
    }
}

impl From<BodePlot> for BodePlotEgui {
    fn from(value: BodePlot) -> Self {
        Self::new(value)
    }
}

impl eframe::App for BodePlotEgui {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.with_layout(Layout::right_to_left(Align::RIGHT), |ui| {
                if ui.button(egui::RichText::new("🏠").size(20.0)).clicked() {
                    *self = self.bode.clone().into();
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let plot_height = ui.available_height() / 2.0;
            ui.vertical(|ui| {
                ui.add_sized(
                    [ui.available_width(), plot_height],
                    |ui: &mut egui::Ui| self.add_mag_plot(ui),
                );

                ui.add_sized(
                    [ui.available_width(), plot_height],
                    |ui: &mut Ui| self.add_phase_plot(ui),
                );
            });
        });
        self.init = true;
    }
}

impl BodePlotEgui {
    fn add_mag_plot(&mut self, ui: &mut Ui) -> egui::Response {
        egui_plot::Plot::new("Magnitude Plot")
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                if !self.init {
                    let bounds = BodePlotEgui::get_mag_bounds(
                        self.bode.options.x_limits,
                        None,
                        self.bode
                            .plot_data
                            .iter()
                            .flat_map(|data| data.mag_phase_freq_points.iter()),
                    );
                    plot_ui.set_plot_bounds(bounds);
                }

                BodePlotEgui::plot_mag_data(plot_ui, &self.bode.plot_data);
                [self.shared_x_limits, self.last_mag_x_limits] =
                    BodePlotEgui::synchronize_x_axes(
                        plot_ui,
                        self.last_mag_x_limits,
                        self.shared_x_limits,
                    );
            })
            .response
    }
    fn plot_mag_data(plot_ui: &mut PlotUi, data: &[BodePlotData]) {
        for data_i in data {
            let plot_points = PlotPoints::from_iter(
                data_i
                    .mag_phase_freq_points
                    .iter()
                    .map(|p| [p[2].log10(), p[0]]),
            );
            let mut line = Line::new(plot_points);
            line = set_name_and_color_line(line, &data_i.name, &data_i.color);
            plot_ui.line(line);
        }
    }

    fn get_mag_bounds<'a, I>(
        x_limits: Option<[f64; 2]>,
        y_limits: Option<[f64; 2]>,
        data: I,
    ) -> PlotBounds
    where
        I: Iterator<Item = &'a [f64; 3]>,
    {
        let xy_iter =
            data.map(|mag_phase_freq| [mag_phase_freq[2], mag_phase_freq[0]]);
        let bounds = get_2d_plot_bounds(x_limits, y_limits, xy_iter);
        let [x_min, y_min] = bounds.min();
        let [x_max, y_max] = bounds.max();
        PlotBounds::from_min_max([x_min.log10(), y_min], [x_max.log10(), y_max])
    }

    fn add_phase_plot(&mut self, ui: &mut Ui) -> egui::Response {
        egui_plot::Plot::new("Phase Plot")
            .legend(Legend::default())
            .show(ui, |plot_ui| {
                if !self.init {
                    let bounds = BodePlotEgui::get_phase_bounds(
                        self.bode.options.x_limits,
                        None,
                        self.bode
                            .plot_data
                            .iter()
                            .flat_map(|data| data.mag_phase_freq_points.iter()),
                    );
                    plot_ui.set_plot_bounds(bounds);
                }
                BodePlotEgui::plot_phase_data(plot_ui, &self.bode.plot_data);
                [self.shared_x_limits, self.last_phase_x_limits] =
                    BodePlotEgui::synchronize_x_axes(
                        plot_ui,
                        self.last_phase_x_limits,
                        self.shared_x_limits,
                    );
            })
            .response
    }

    fn get_phase_bounds<'a, I>(
        x_limits: Option<[f64; 2]>,
        y_limits: Option<[f64; 2]>,
        data: I,
    ) -> PlotBounds
    where
        I: Iterator<Item = &'a [f64; 3]>,
    {
        let xy_iter =
            data.map(|mag_phase_freq| [mag_phase_freq[2], mag_phase_freq[1]]);
        let bounds = get_2d_plot_bounds(x_limits, y_limits, xy_iter);
        let [x_min, y_min] = bounds.min();
        let [x_max, y_max] = bounds.max();
        PlotBounds::from_min_max([x_min.log10(), y_min], [x_max.log10(), y_max])
    }

    fn plot_phase_data(plot_ui: &mut PlotUi, data: &[BodePlotData]) {
        for data_i in data {
            let plot_points = PlotPoints::from_iter(
                data_i
                    .mag_phase_freq_points
                    .iter()
                    .map(|p| [p[2].log10(), p[1]]),
            );
            let mut line: Line<'_> = Line::new(plot_points);
            line = set_name_and_color_line(line, &data_i.name, &data_i.color);
            plot_ui.line(line);
        }
    }
    fn synchronize_x_axes(
        plot_ui: &mut PlotUi,
        last_x_limits: [f64; 2],
        shared_x_limits: [f64; 2],
    ) -> [[f64; 2]; 2] {
        let current_x_limits = bounds_to_x_limits(&plot_ui.plot_bounds());
        let mut new_shared_x_limits = shared_x_limits;
        if last_x_limits != current_x_limits {
            new_shared_x_limits = current_x_limits;
        } else if current_x_limits != shared_x_limits {
            let mut bounds = plot_ui.plot_bounds();
            bounds = bounds_set_x_limits(&bounds, shared_x_limits);
            plot_ui.set_plot_bounds(bounds);
        }
        [new_shared_x_limits, current_x_limits]
    }
}

impl Plot for BodePlot {
    fn show(
        self,
        width: usize,
        height: usize,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        BodePlotEgui::new(self).show(width, height, name)
    }
}

fn get_2d_plot_bounds<I>(
    x_limits: Option<[f64; 2]>,
    y_limits: Option<[f64; 2]>,
    xy_iter: I,
) -> PlotBounds
where
    I: Iterator<Item = [f64; 2]>,
{
    let [x_min_init, x_max_init] =
        x_limits.unwrap_or([-f64::INFINITY, f64::INFINITY]);
    let [y_min_init, y_max_init] =
        y_limits.unwrap_or([-f64::INFINITY, f64::INFINITY]);
    let filtered_data = xy_iter.filter(|c| {
        c[0] >= x_min_init
            && c[0] <= x_max_init
            && c[1] >= y_min_init
            && c[1] <= y_max_init
    });

    let [mut x_min, mut x_max] = [f64::INFINITY, -f64::INFINITY];
    let [mut y_min, mut y_max] = [f64::INFINITY, -f64::INFINITY];

    for c in filtered_data {
        x_min = x_min.min(c[0]);
        x_max = x_max.max(c[0]);
        y_min = y_min.min(c[1]);
        y_max = y_max.max(c[1]);
    }

    [x_min, x_max] = x_limits.unwrap_or([x_min, x_max]);
    [y_min, y_max] = y_limits.unwrap_or([y_min, y_max]);

    PlotBounds::from_min_max([x_min, y_min], [x_max, y_max])
}

fn set_name_and_color_line<'a>(
    mut line: Line<'a>,
    name: &str,
    color: &Option<RGBAColor>,
) -> Line<'a> {
    if !name.is_empty() {
        line = line.name(name);
    }
    if let Some(color) = color {
        line = line.color(color.to_egui());
    }
    line
}

fn bounds_to_x_limits(bounds: &PlotBounds) -> [f64; 2] {
    [bounds.min()[0], bounds.max()[0]]
}

fn bounds_set_x_limits(bounds: &PlotBounds, x_limits: [f64; 2]) -> PlotBounds {
    PlotBounds::from_min_max([x_limits[0], bounds.min()[1]], [
        x_limits[1],
        bounds.max()[1],
    ])
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
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.with_layout(Layout::right_to_left(Align::RIGHT), |ui| {
                if ui.button(egui::RichText::new("🏠").size(20.0)).clicked() {
                    *self = self.nyq.clone().into();
                }
            });
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            let plot =
                egui_plot::Plot::new("Nyquist Plot").legend(Legend::default());
            plot.show(ui, |plot_ui| {
                if !self.init {
                    let data_iter = self
                        .nyq
                        .plot_data
                        .iter()
                        .flat_map(|nyq_data| nyq_data.data_points.iter());

                    let bounds_init = NyquistPlotEgui::get_plot_bounds(
                        self.nyq.options.x_limits,
                        self.nyq.options.y_limits,
                        data_iter,
                    );
                    plot_ui.set_plot_bounds(bounds_init);
                }
                NyquistPlotEgui::plot_data(plot_ui, &self.nyq.plot_data);
            });
        });
        self.init = true;
    }
}

impl NyquistPlotEgui {
    fn get_plot_bounds<'a, I>(
        x_limits: Option<[f64; 2]>,
        y_limits: Option<[f64; 2]>,
        data: I,
    ) -> PlotBounds
    where
        I: Iterator<Item = &'a Complex64>,
    {
        let xy_iter = data.map(|c| [c.re, c.im]);
        get_2d_plot_bounds(x_limits, y_limits, xy_iter)
    }

    fn plot_data(plot_ui: &mut PlotUi, data: &[NyquistPlotData]) {
        for data_i in data {
            let plot_points = PlotPoints::from_iter(
                data_i.data_points.iter().map(|c| [c.re, c.im]),
            );
            let mut line = Line::new(plot_points);
            line = set_name_and_color_line(line, &data_i.name, &data_i.color);
            plot_ui.line(line);
        }
    }
}

impl Plot for NyquistPlot {
    fn show(
        self,
        width: usize,
        height: usize,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        NyquistPlotEgui::new(self).show(width, height, name)
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
    use egui_kittest::Harness;

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

        let app = BodePlotEgui::new(bodeplot);
        let mut harness = Harness::new_eframe(|_| app);
        harness.run_ok().unwrap();
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

        let mut harness =
            Harness::new_eframe(|_| NyquistPlotEgui::new(nyq_plot));
        harness.run_ok().unwrap();
    }
}
