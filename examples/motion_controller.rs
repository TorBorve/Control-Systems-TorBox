use control_systems_torbox::{
    BodePlot, BodePlotOptions, Continuous, FrequencyResponse, NyquistPlot,
    NyquistPlotOptions, Plot, Ss, Tf, analysis::frequency_response::lin_space,
};

fn main() {
    // This example shows how we can create a plant, design a controller for the
    // plant and analyse the resulting closed loop system.

    // Plant is integrator, 1/s, and a power converter with a bandwith of 10
    // rad/s, and a time delay of 0.5s approximated with a second order pade
    // approximation
    let s = Tf::s();
    let plant = 1.0 / &s * 50.0 / (&s + 50.0) * Tf::pade(0.12, 3);
    // let plant = plant.to_ss().unwrap();

    // We will use a PI controller: Kp*(1 + 1/T_i*1/s)
    // We set the integration time to 1 second and find the largest Kp such that
    // Ms = ||S(s)||_inf = peak of sensitity function is less than two.
    let ms_limit = 2.0;
    let t_i = 10.0;
    let kps = lin_space(100.0, 0.001, 100); // Kp to test
    let mut controller = None;
    for kp in kps {
        let controller_i = kp * (1.0 + 1.0 / t_i * 1.0 / &s); //.to_ss().unwrap();
        let open_loop = controller_i.series(&plant);
        // S = 1/(1 + L)
        let sensitivity_func =
            Ss::new_from_scalar(1.0).feedback(&open_loop.to_ss().unwrap());
        if !sensitivity_func.is_stable() {
            continue;
        }

        let max_sensitivity = sensitivity_func.norm_hinf().unwrap();

        // println!("Ms = {}", max_sensitivity);
        // println!("Stable: {}", sensitivity_func.is_stable());
        let mut bode_plot = BodePlot::new(BodePlotOptions::default());
        bode_plot.add_system(sensitivity_func.bode(0.1, 100.).into());
        bode_plot.add_system(
            Tf::<f64, Continuous>::new_from_scalar(ms_limit)
                .bode(0.1, 100.)
                .into(),
        );
        // uncomment if you want to see the progression of the search
        // bode_plot.show(600, 400, "Candidate Sensitivity Function").unwrap();
        if max_sensitivity < ms_limit {
            println!("Found Kp = {}, which gives Ms = {}", kp, max_sensitivity);
            controller = Some(controller_i);
            break;
        }
    }
    let controller = controller.unwrap();

    // Lets look at the properties of the controller and closed loop system

    println!("Open Loop transfer function");
    let open_loop = controller.series(&plant);
    let mut bode_plot = BodePlot::new(BodePlotOptions::default());
    bode_plot.add_system(open_loop.bode(0.1, 100.0).into());
    bode_plot.show(600, 400, "Open Loop").unwrap();

    println!("Sensitivity Function bode plot");
    let sensitivity_func = Tf::new_from_scalar(1.0).feedback(&open_loop);
    bode_plot = BodePlot::new(BodePlotOptions::default());
    bode_plot.add_system(sensitivity_func.bode(0.1, 100.).into());
    bode_plot.add_system(
        Tf::<f64, Continuous>::new_from_scalar(ms_limit)
            .bode(0.1, 100.)
            .into(),
    );
    bode_plot.show(600, 400, "Sensitivity Function").unwrap();

    println!("Nyquist plot of open loop");
    let mut nyq_plot = NyquistPlot::new(NyquistPlotOptions::default());
    nyq_plot.add_system(open_loop.nyquist(0.1, 100.).into());
    nyq_plot.show(400, 400, "Nyquist Open Loop").unwrap();

    println!("Tracking Function bode plot");
    let tracking_func = open_loop.feedback(&Tf::new_from_scalar(1.0));
    let mut bode_plot = BodePlot::new(BodePlotOptions::default());
    bode_plot.add_system(tracking_func.bode(0.1, 100.).into());
    bode_plot
        .show(600, 400, "Tracking Function Bode Plot")
        .unwrap();

    println!("Location of poles: [re, im]");
    let poles = tracking_func.to_ss().unwrap().poles();
    for pole in poles {
        print!("[{:.2}, {:.2}] ", pole.re, pole.im);
    }
    println!();
}
