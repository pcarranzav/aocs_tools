//! This example runs a simulation of a simple detumble controller using Modified Rodrigues Parameters
use aocs_tools::{control::{mrp_zero_torque, MRPDetumbleController, MRPController}, math::Vector6, plotting::plot_mrp_simulation};

fn main() {
    let inertia = [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 2.0]];
    let initial_mrp = [0.0, 0.0, 0.0];
    let initial_omega = [1.0, 0.5, 0.1];
    let controller = MRPDetumbleController::new(0.5);
    let control = |time: f64, state: &Vector6| controller.control_torque(time, state);
    let (x, t) = aocs_tools::simulation::mrp_attitude_simulator(
        0.001,
        40.0,
        inertia,
        initial_mrp,
        initial_omega,
        control,
        mrp_zero_torque,
    );
    plot_mrp_simulation(&x, &t);
}
