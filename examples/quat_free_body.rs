//! This example runs a simulation of a free body in 3D space, using a quaternion attitude representation
use aocs_tools::{control::quat_zero_torque, plotting::plot_quat_simulation};
use quaternion::Quaternion;

fn main() {
    let inertia = [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 2.0]];
    let initial_q: Quaternion<f64> = (1.0, [0.0, 0.0, 0.0]);
    let initial_omega = [1.0, 0.5, 0.1];
    let (x, t) = aocs_tools::simulation::quat_attitude_simulator(
        0.001,
        10.0,
        inertia,
        initial_q,
        initial_omega,
        quat_zero_torque,
        quat_zero_torque,
    );
    plot_quat_simulation(&x, &t);
}
