//! Simulation module (WIP)

use quaternion::Quaternion;
use vecmath::{Matrix3, Vector3};

use crate::{
    kinematics::{self, mrp_differential_equation_matrix},
    math::{self, Vector6, VectorN, Vector7}, kinetics,
};

/// Generic RK4 simulator.
/// It takes in the state derivative function, the initial conditions, the step size, the final time, and an optional postprocessing function
/// The postprocessing function sanitizes the state as in xn = p(xn); e.g. constraining its module or normalizing it.
/// The function returns the state and time vectors.
/// Note the state can be any size of vector, as it depends on the attitude representation or whatever we're trying to simulate,
/// the vector type must implement math::VectorN (addition and multiply by scalar).
pub fn rk4_simulator<T, S, P>(
    f: S,
    x0: T,
    step: f64,
    duration: f64,
    p: Option<P>,
) -> (Vec<T>, Vec<f64>)
where
    T: VectorN,
    S: Fn(f64, &T) -> T,
    P: Fn(&T) -> T,
{
    let mut t = 0.0;
    let mut x = vec![x0];
    let mut t_vec = vec![t];
    while t < duration {
        let x_last = x.last().unwrap(); // Must exist because we're always including an element
        let k1 = f(t, &x_last).vec_scale(step);
        let t_mid = t + step / 2.0;
        let k2 = f(t_mid, &x_last.vec_add(&k1.vec_scale(0.5))).vec_scale(step);
        let k3 = f(t_mid, &x_last.vec_add(&k2.vec_scale(0.5))).vec_scale(step);
        t += step;
        let k4 = f(t, &x_last.vec_add(&k3)).vec_scale(step);

        let k2 = k2.vec_scale(2.0);
        let k3 = k3.vec_scale(2.0);

        let k_total = k1
            .vec_add(&k2)
            .vec_add(&k3)
            .vec_add(&k4)
            .vec_scale(1.0 / 6.0);
        let mut xn = x_last.vec_add(&k_total);
        if let Some(p) = &p {
            xn = p(&xn);
        }
        x.push(xn);
        t_vec.push(t);
    }
    (x, t_vec)
}

/// Simulator for a rigid spacecraft's attitude using the MRP representation
/// The state representation is a 6D vector, with the first 3 elements being the MRP attitude and the last 3 the angular velocity in the B frame
pub fn mrp_attitude_simulator<U, D>(
    step: f64,
    duration: f64,
    inertia: Matrix3<f64>,
    initial_mrp: Vector3<f64>,
    initial_omega: Vector3<f64>,
    control_torque: U,
    disturbance_torque: D,
) -> (Vec<Vector6>, Vec<f64>)
where
    U: Fn(f64, &Vector6) -> Vector3<f64>,
    D: Fn(f64, &Vector6) -> Vector3<f64>,
{
    let x0 = [
        initial_mrp[0],
        initial_mrp[1],
        initial_mrp[2],
        initial_omega[0],
        initial_omega[1],
        initial_omega[2],
    ];
    let inertia_inverse = vecmath::mat3_inv(inertia);

    let x_dot = |t: f64, x: &Vector6| {
        let mrp = [x[0], x[1], x[2]];
        let omega = [x[3], x[4], x[5]];
        let mrp_dot = vecmath::row_mat3_transform(
            math::mat3_scale(mrp_differential_equation_matrix(mrp), 0.25),
            omega,
        );
        let torque = vecmath::vec3_add(control_torque(t, x), disturbance_torque(t, x));
        let omega_dot = kinetics::rigid_body_omega_dot(inertia, omega, torque, Some(inertia_inverse));
        [
            mrp_dot[0],
            mrp_dot[1],
            mrp_dot[2],
            omega_dot[0],
            omega_dot[1],
            omega_dot[2],
        ]
    };
    let p = |x: &Vector6| {
        let mrp = [x[0], x[1], x[2]];
        let mrp_norm_squared = vecmath::vec3_square_len(mrp);
        if mrp_norm_squared > 1.0 {
            let mrp_normalized = kinematics::mrp_shadow(mrp);
            [
                mrp_normalized[0],
                mrp_normalized[1],
                mrp_normalized[2],
                x[3],
                x[4],
                x[5],
            ]
        } else {
            *x
        }
    };
    rk4_simulator(x_dot, x0, step, duration, Some(p))
}

/// Simulator for a rigid spacecraft's attitude using a quaternion representation
/// The state representation is a 7D vector, with the first 4 elements being the quaternion attitude and the last 3 the angular velocity in the B frame
pub fn quat_attitude_simulator<U, D>(
    step: f64,
    duration: f64,
    inertia: Matrix3<f64>,
    initial_quat: Quaternion<f64>,
    initial_omega: Vector3<f64>,
    control_torque: U,
    disturbance_torque: D,
) -> (Vec<Vector7>, Vec<f64>)
where
    U: Fn(f64, &Vector7) -> Vector3<f64>,
    D: Fn(f64, &Vector7) -> Vector3<f64>,
{
    let x0 = [
        initial_quat.0,
        initial_quat.1[0],
        initial_quat.1[1],
        initial_quat.1[2],
        initial_omega[0],
        initial_omega[1],
        initial_omega[2],
    ];
    let inertia_inverse = vecmath::mat3_inv(inertia);

    let x_dot = |t: f64, x: &Vector7| {
        let q: Quaternion<f64> = (x[0], [x[1], x[2], x[3]]);
        let omega = [x[4], x[5], x[6]];
        let q_dot = kinematics::omega_to_qdot(q, omega);
        let torque = vecmath::vec3_add(control_torque(t, x), disturbance_torque(t, x));
        let omega_dot = kinetics::rigid_body_omega_dot(inertia, omega, torque, Some(inertia_inverse));
        [
            q_dot.0,
            q_dot.1[0],
            q_dot.1[1],
            q_dot.1[2],
            omega_dot[0],
            omega_dot[1],
            omega_dot[2],
        ]
    };
    let p = |x: &Vector7| {
        let q: Quaternion<f64> = (x[0], [x[1], x[2], x[3]]);
        let q_normalized = math::quat_normalize(q);
        [
            q_normalized.0,
            q_normalized.1[0],
            q_normalized.1[1],
            q_normalized.1[2],
            x[4],
            x[5],
            x[6],
        ]
    };
    rk4_simulator(x_dot, x0, step, duration, Some(p))
}
