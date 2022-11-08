//! Kinetics module (WIP)

use vecmath::{Matrix3, Vector3};

use crate::math;

/// Calculate the time derivative of the angular velocity based on inertia and torque
/// Note the function accepts an optional inverse for the inertia tensor, to avoid recomputing it
pub fn rigid_body_omega_dot(
    inertia: Matrix3<f64>,
    omega: Vector3<f64>,
    torque: Vector3<f64>,
    inertia_inverse: Option<Matrix3<f64>>,
) -> Vector3<f64> {
    let inertia_inverse = inertia_inverse.unwrap_or_else(|| vecmath::mat3_inv(inertia));
    let omega_tilde = math::vec3_tilde(omega);
    let omega_tilde_inertia_omega =
        vecmath::row_mat3_transform(omega_tilde, vecmath::row_mat3_transform(inertia, omega));
    let inertia_omega_dot = vecmath::vec3_sub(
        torque,
        omega_tilde_inertia_omega,
    );
    vecmath::row_mat3_transform(inertia_inverse, inertia_omega_dot)
}
