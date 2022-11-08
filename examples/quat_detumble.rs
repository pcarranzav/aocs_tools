//! This example runs a simulation of a simple detumble controller using quaternions
use aocs_tools::{control::{quat_zero_torque, QuatDetumbleController, QuatController}, math::{Vector7}, plotting::plot_quat_simulation};
use quaternion::Quaternion;

fn main() {
    let inertia = [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 2.0]];
    let initial_q: Quaternion<f64> = (1.0, [0.0, 0.0, 0.0]);
    let initial_omega = [1.0, 0.5, 0.1];
    let controller = QuatDetumbleController::new(0.5);
    let control = |time: f64, state: &Vector7| controller.control_torque(time, state);
    let (x, t) = aocs_tools::simulation::quat_attitude_simulator(
        0.001,
        40.0,
        inertia,
        initial_q,
        initial_omega,
        control,
        quat_zero_torque,
    );
    plot_quat_simulation(&x, &t);
}
