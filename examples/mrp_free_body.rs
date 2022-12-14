//! This example runs a simulation of a free body in 3D space, using a Modified Rodrigues Parameters attitude representation
use aocs_tools::{control::mrp_zero_torque, plotting::plot_mrp_simulation};

fn main() {
    let inertia = [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 2.0]];
    let initial_mrp = [0.0, 0.0, 0.0];
    let initial_omega = [1.0, 0.5, 0.1];
    let (x, t) = aocs_tools::simulation::mrp_attitude_simulator(
        0.001,
        10.0,
        inertia,
        initial_mrp,
        initial_omega,
        mrp_zero_torque,
        mrp_zero_torque,
    );
    plot_mrp_simulation(&x, &t);
}
