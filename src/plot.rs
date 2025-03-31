use core::f64;

use egui::{Align, Color32, Layout, Ui};
use egui_plot::{Legend, Line, PlotBounds, PlotPoints, PlotUi};
use num_complex::Complex64;

/// A struct representing a color in the RGBA color space.
///
/// This struct holds four components: red (`r`), green (`g`), blue (`b`), and
/// alpha (`a`) channels, where each channel is an 8-bit unsigned integer,
/// meaning it ranges from 0 to 255.
///
/// # Examples
///
/// ```rust
/// use control_systems_torbox::RGBAColor;
/// let red = RGBAColor::RED;
/// assert_eq!(red.r, 255);
/// assert_eq!(red.g, 0);
/// assert_eq!(red.b, 0);
/// assert_eq!(red.a, 255);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RGBAColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl RGBAColor {
    /// Common colors represented by RGBAColor constants.
    ///
    /// These are predefined color constants for commonly used colors.
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
    /// Creates a new `RGBAColor` instance with the specified red, green, blue,
    /// and alpha values.
    ///
    /// # Parameters
    /// - `red`: The red channel value (0-255).
    /// - `green`: The green channel value (0-255).
    /// - `blue`: The blue channel value (0-255).
    /// - `alpha`: The alpha (transparency) channel value (0-255).
    pub fn new(red: u8, green: u8, blue: u8, alpha: u8) -> Self {
        Self {
            r: red,
            g: green,
            b: blue,
            a: alpha,
        }
    }

    /// Converts the `RGBAColor` to a `Color32` type used in the `egui` library.
    ///
    /// # Returns
    /// A `Color32` representation of the color.
    fn to_egui(self) -> Color32 {
        Color32::from_rgba_premultiplied(self.r, self.g, self.b, self.a)
    }
}

/// The `Plot` trait defines a common interface for types that can display a
/// plot. Types implementing this trait should be able to show a plot in a
/// graphical window.
pub trait Plot {
    /// Displays the plot in a new window with the given dimensions and name.
    ///
    /// # Arguments
    ///
    /// * `width` - The width of the plot window in pixels.
    /// * `height` - The height of the plot window in pixels.
    /// * `name` - The name of the plot window.
    ///
    /// # Returns
    ///
    /// Returns a `Result` indicating success or failure.
    fn show(
        self,
        width: usize,
        height: usize,
        name: &str,
    ) -> Result<(), Box<dyn std::error::Error + 'static>>;
}

impl<T: eframe::App> Plot for T {
    /// Implements the `Plot` trait for types that implement the `eframe::App`
    /// trait.
    ///
    /// This function opens a new window with the given dimensions and runs the
    /// application that implements the `eframe::App` trait.
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

/// `BodePlotData` represents the data for a Bode plot, including magnitude,
/// phase, and frequency points.
#[derive(Clone)]
pub struct BodePlotData {
    mag_phase_freq_points: Vec<[f64; 3]>,
    name: String,
    color: Option<RGBAColor>,
}

impl BodePlotData {
    /// Creates a new `BodePlotData` with the given magnitude, phase, and
    /// frequency points.
    ///
    /// # Arguments
    ///
    /// * `mag_phase_freq_points` - A vector of 3-element arrays representing
    ///   the magnitude, phase, and frequency points for the Bode plot.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotData`.
    pub fn new(mag_phase_freq_points: Vec<[f64; 3]>) -> Self {
        Self {
            mag_phase_freq_points,
            name: "".to_string(),
            color: None,
        }
    }

    /// Sets the name for the Bode plot data.
    ///
    /// # Arguments
    ///
    /// * `name` - The name to set for the Bode plot data.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotData` with the updated name.
    pub fn set_name(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }

    /// Sets the color for the Bode plot data.
    ///
    /// # Arguments
    ///
    /// * `color` - The color to set for the Bode plot data.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotData` with the updated color.
    pub fn set_color(mut self, color: RGBAColor) -> Self {
        self.color = Some(color);
        self
    }
}

impl From<Vec<[f64; 3]>> for BodePlotData {
    /// Converts a vector of 3-element arrays to a `BodePlotData` instance.
    ///
    /// # Arguments
    ///
    /// * `value` - A vector of 3-element arrays representing magnitude, phase,
    ///   and frequency points.
    ///
    /// # Returns
    ///
    /// A `BodePlotData` instance.
    fn from(value: Vec<[f64; 3]>) -> Self {
        BodePlotData::new(value)
    }
}

/// `BodePlotOptions` holds the options for configuring a Bode plot, such as the
/// title and axis limits.
#[derive(Default, Clone)]
pub struct BodePlotOptions {
    title: String,
    x_limits: Option<[f64; 2]>,
}

impl BodePlotOptions {
    /// Creates a new `BodePlotOptions` instance with the given title.
    ///
    /// # Arguments
    ///
    /// * `title` - The title to set for the Bode plot.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotOptions`.
    pub fn new(title: &str) -> Self {
        Self {
            title: title.to_string(),
            x_limits: None,
        }
    }

    /// Sets the title for the Bode plot options.
    ///
    /// # Arguments
    ///
    /// * `title` - The title to set for the Bode plot.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotOptions` with the updated title.
    pub fn set_title(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }

    /// Sets the x-axis limits for the Bode plot options.
    ///
    /// # Arguments
    ///
    /// * `x_limits` - The x-axis limits as a 2-element array `[min, max]`.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotOptions` with the updated x-axis limits.
    pub fn set_x_limits(mut self, x_limits: [f64; 2]) -> Self {
        self.x_limits = Some(x_limits);
        self
    }
}

/// `BodePlot` holds the configuration options and data for a Bode plot.
#[derive(Default, Clone)]
pub struct BodePlot {
    options: BodePlotOptions,
    plot_data: Vec<BodePlotData>,
}

impl BodePlot {
    /// Creates a new `BodePlot` with the given options.
    ///
    /// # Arguments
    ///
    /// * `options` - The options to configure the Bode plot.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlot`.
    pub fn new(options: BodePlotOptions) -> Self {
        Self {
            options,
            plot_data: vec![],
        }
    }

    // Adds a system's Bode plot data to the `BodePlot`.
    ///
    /// # Arguments
    ///
    /// * `data` - The `BodePlotData` for the system.
    ///
    /// # Returns
    ///
    /// A mutable reference to the `BodePlot` instance for method chaining.
    pub fn add_system(&mut self, data: BodePlotData) -> &mut Self {
        self.plot_data.push(data);
        self
    }
}

/// `BodePlotEgui` is a wrapper for displaying the `BodePlot` using the `eframe`
/// GUI framework.
#[derive(Clone)]
pub struct BodePlotEgui {
    bode: BodePlot,
    init: bool,
    shared_x_limits: [f64; 2],
    last_mag_x_limits: [f64; 2],
    last_phase_x_limits: [f64; 2],
}

impl BodePlotEgui {
    /// Creates a new `BodePlotEgui` instance with the given `BodePlot` data.
    ///
    /// # Arguments
    ///
    /// * `bode` - The `BodePlot` instance to be wrapped by the `BodePlotEgui`.
    ///
    /// # Returns
    ///
    /// A new instance of `BodePlotEgui`.
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
    /// Converts a `BodePlot` instance into a `BodePlotEgui` instance.
    fn from(value: BodePlot) -> Self {
        Self::new(value)
    }
}

impl eframe::App for BodePlotEgui {
    /// Updates the `BodePlotEgui` GUI on each frame.
    ///
    /// This function is called to update the user interface and plot the
    /// magnitude and phase plots.
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

/// Implements the `Plot` trait for the `BodePlot` struct.
///
/// This implementation defines how to display a `BodePlot` using the `Plot`
/// trait's `show` method, which delegates the rendering of the plot to the
/// `BodePlotEgui` struct for use with the `eframe` framework.
///
/// # Arguments
///
/// * `self` - The `BodePlot` instance that is being displayed.
/// * `width` - The width of the plot window in pixels.
/// * `height` - The height of the plot window in pixels.
/// * `name` - A string representing the name of the plot window, used for
///   identification in the UI.
///
/// # Returns
///
/// This method returns a `Result<(), Box<dyn std::error::Error + 'static>>`,
/// which represents the success or failure of showing the plot. If successful,
/// the result is `Ok(())`, and if there is an error, it returns a boxed error
/// type.
///
/// # Example
///
/// ```rust,no_run
/// use control_systems_torbox::*;
/// let bode_data = vec![
///     [1.0, 2.0, 3.0],
///     [1.1, 2.1, 3.1],
///     [1.2, 2.2, 3.2],
/// ];
///
/// let mut bode = BodePlot::new(BodePlotOptions::new("Sample Bode Plot"));
/// bode.add_system(BodePlotData::new(bode_data));
///
/// bode.show(800, 600, "Bode Plot Example").unwrap();
/// ```
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
    PlotBounds::from_min_max(
        [x_limits[0], bounds.min()[1]],
        [x_limits[1], bounds.max()[1]],
    )
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Nyquist
//////////////////////////////////////////////////////////////////////////////////////////////

/// A struct representing the data points for a Nyquist plot.
///
/// This struct stores the data points for the Nyquist plot, as well as optional
/// metadata such as the name of the system and the color for plotting.
///
/// # Fields
///
/// - `data_points`: A `Vec<Complex64>` that contains the data points to be
///   plotted.
/// - `name`: A `String` representing the name of the system or dataset. It is
///   used for labeling the plot.
/// - `color`: An optional `RGBAColor` that defines the color of the plot line.
#[derive(Clone, Debug)]
pub struct NyquistPlotData {
    data_points: Vec<Complex64>,
    name: String,
    color: Option<RGBAColor>,
}

impl NyquistPlotData {
    /// Creates a new `NyquistPlotData` instance with the given data points.
    ///
    /// # Arguments
    /// - `data_points`: A vector of complex numbers representing the Nyquist
    ///   plot data points.
    ///
    /// # Returns
    /// A new instance of `NyquistPlotData`.
    pub fn new(data_points: Vec<Complex64>) -> Self {
        Self {
            data_points,
            name: "".to_string(),
            color: None,
        }
    }

    /// Sets the name of the Nyquist plot data.
    ///
    /// # Arguments
    /// - `name`: The name to be set for the Nyquist plot data.
    ///
    /// # Returns
    /// A new instance of `NyquistPlotData` with the name set.
    pub fn set_name(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }

    /// Sets the color for the Nyquist plot data.
    ///
    /// # Arguments
    /// - `color`: The `RGBAColor` to be used for plotting.
    ///
    /// # Returns
    /// A new instance of `NyquistPlotData` with the color set.
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

/// A struct representing the options for a Nyquist plot.
///
/// This struct holds the configuration settings for the Nyquist plot, such as
/// title, axis limits, etc.
///
/// # Fields
///
/// - `title`: The title of the plot.
/// - `x_limits`: Optional X-axis limits as an array `[min, max]`.
/// - `y_limits`: Optional Y-axis limits as an array `[min, max]`.
#[derive(Default, Clone, Debug)]
pub struct NyquistPlotOptions {
    title: String,
    x_limits: Option<[f64; 2]>,
    y_limits: Option<[f64; 2]>,
}

impl NyquistPlotOptions {
    /// Sets the title for Nyquist plot.
    pub fn set_title(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }

    /// Sets the x-limits for Nyquist plot.
    pub fn set_x_limits(mut self, x_limits: [f64; 2]) -> Self {
        self.x_limits = Some(x_limits);
        self
    }

    /// Sets the y-limits for Nyquist plot.
    pub fn set_y_limits(mut self, y_limits: [f64; 2]) -> Self {
        self.y_limits = Some(y_limits);
        self
    }
}

/// A struct representing the Nyquist plot.
///
/// This struct contains the data and options to plot a Nyquist plot, including
/// the `NyquistPlotOptions` and the actual data in `NyquistPlotData` format.
///
/// # Fields
/// - `options`: The configuration options for the Nyquist plot (title, axis
///   limits).
/// - `plot_data`: A vector of `NyquistPlotData` representing the systems or
///   datasets to plot.
#[derive(Debug, Default, Clone)]
pub struct NyquistPlot {
    options: NyquistPlotOptions,
    plot_data: Vec<NyquistPlotData>,
}

impl NyquistPlot {
    /// Creates a new `NyquistPlot` with the given options.
    ///
    /// # Arguments
    /// - `options`: The configuration options for the plot.
    ///
    /// # Returns
    /// A new instance of `NyquistPlot`.
    pub fn new(options: NyquistPlotOptions) -> Self {
        Self {
            options,
            plot_data: vec![],
        }
    }

    // Adds a system's data to the Nyquist plot.
    ///
    /// # Arguments
    /// - `data`: The data for a system to be plotted.
    ///
    /// # Returns
    /// A mutable reference to the `NyquistPlot` instance for method chaining.
    pub fn add_system(&mut self, data: NyquistPlotData) -> &mut Self {
        self.plot_data.push(data);
        self
    }
}

/// A helper struct for rendering the Nyquist plot using the `eframe` framework.
///
/// This struct wraps a `NyquistPlot` and handles the display logic using the
/// `eframe::App` trait.
///
/// # Fields
/// - `nyq`: The `NyquistPlot` instance to be displayed.
/// - `init`: A flag to ensure that the plot bounds are set only once.
pub struct NyquistPlotEgui {
    nyq: NyquistPlot,
    init: bool,
}

impl NyquistPlotEgui {
    /// Creates a new `NyquistPlotEgui` instance.
    ///
    /// # Arguments
    /// - `nyq`: The `NyquistPlot` instance to be rendered.
    ///
    /// # Returns
    /// A new instance of `NyquistPlotEgui`.
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
    /// Implements the `Plot` trait for `NyquistPlot`.
    ///
    /// This implementation allows the `NyquistPlot` to be displayed in a window
    /// by utilizing the `NyquistPlotEgui` struct.
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
        analysis::frequency_response::{bode, nyquist},
        systems::Tf,
        utils::traits::Continuous,
    };
    use egui_kittest::Harness;

    #[test]
    fn colors() {
        let red = RGBAColor::RED;
        let red_new = RGBAColor::new(255, 0, 0, 255);
        assert_eq!(red_new, red);
    }

    #[test]
    fn bodeplot() {
        let bodeopts = BodePlotOptions::new("ERROR WRONG TITLE")
            .set_title("Bodeplot title")
            .set_x_limits([0.1, 100.]);
        let mut bodeplot = BodePlot::new(bodeopts);

        let sys: Tf<f64, Continuous> = Tf::new(&[0.0, 1.0], &[1.0, 1.0]);
        // let (mag, phase, freq) = bode(sys, 0.1, 100.);
        // let mag_phase_freq_vec: Vec<[f64; 3]> = mag
        //     .iter()
        //     .zip(phase.iter().zip(freq.iter()))
        //     .map(|(&mag, (&phase, &freq))| [mag, phase, freq])
        //     .collect();
        let mag_phase_freq_vec = bode(sys, 0.1, 100.);
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
        let mag_phase_freq_vec = bode(sys2, 0.1, 100.);
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

        let mut nyq_opts =
            NyquistPlotOptions::default().set_title("Nyquist TITLE");
        nyq_opts = nyq_opts.set_x_limits([-5., 1.]).set_y_limits([-10., 10.]);
        let mut nyq_plot = NyquistPlot::new(nyq_opts);
        let nyq_data = nyquist(sys, 0.01, 100.);
        nyq_plot.add_system(nyq_data.into());

        let data = NyquistPlotData::new(nyquist(sys_p1, 0.01, 100.))
            .set_name("Pade 1")
            .set_color(RGBAColor::ORANGE);
        nyq_plot.add_system(data);
        let data = NyquistPlotData::new(nyquist(sys_p2, 0.01, 100.))
            .set_name("Pade 2");
        nyq_plot.add_system(data);

        let mut harness =
            Harness::new_eframe(|_| NyquistPlotEgui::new(nyq_plot));
        harness.run_ok().unwrap();
    }
}
