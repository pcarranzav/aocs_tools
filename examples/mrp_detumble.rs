use aocs_tools::{control::zero_torque, math::Vector6, plotting::plot_mrp_simulation};
use vecmath::Vector3;

// TODO: move the detumble function to the control module
// I need to figure out a nice way to parameterize the control functions.
// Maybe create controller struct types?

fn detumble(_: f64, x: &Vector6) -> Vector3<f64> {
    let gain = 0.6;
    [-gain * x[3], -gain * x[4], -gain * x[5]]
}

fn main() {
    let inertia = [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 2.0]];
    let initial_mrp = [0.0, 0.0, 0.0];
    let initial_omega = [1.0, 0.5, 0.1];
    let (x, t) = aocs_tools::simulation::mrp_attitude_simulator(
        0.001,
        40.0,
        inertia,
        initial_mrp,
        initial_omega,
        detumble,
        zero_torque,
    );
    plot_mrp_simulation(&x, &t);
}
