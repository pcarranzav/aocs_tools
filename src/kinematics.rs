//! Kinematics module (WIP)
//! Provides tools to describe rotations/attitudes and to
//! convert between different representations of rotations.

use quaternion::Quaternion;
use vecmath::{Vector3, Matrix3, Matrix4};

use crate::math;

/// Compute the matrix for a rotation around the third axis (z)
pub fn euler3_rotation_degrees(z_angle: f64) -> Matrix3<f64> {
    let c = z_angle.to_radians().cos();
    let s = z_angle.to_radians().sin();
    [
        [c, s, 0.0],
        [-s, c, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// Compute the matrix for a rotation around the second axis (y)
pub fn euler2_rotation_degrees(y_angle: f64) -> Matrix3<f64> {
    let c = y_angle.to_radians().cos();
    let s = y_angle.to_radians().sin();
    [
        [c, 0.0, -s],
        [0.0, 1.0, 0.0],
        [s, 0.0, c],
    ]
}

/// Compute the matrix for a rotation around the first axis (x)
pub fn euler1_rotation_degrees(x_angle: f64) -> Matrix3<f64> {
    let c = x_angle.to_radians().cos();
    let s = x_angle.to_radians().sin();
    [
        [1.0, 0.0, 0.0],
        [0.0, c, s],
        [0.0, -s, c],
    ]
}

/// Compute the rotation matrix for a set of 3-2-1 Euler angles (yaw, pitch, roll)
pub fn euler321_to_dcm(theta: &Vector3<f64>) -> Matrix3<f64> {
    let r3 = euler3_rotation_degrees(theta[0]);
    let r2 = euler2_rotation_degrees(theta[1]);
    let r1 = euler1_rotation_degrees(theta[2]);

    let r1r2 = vecmath::row_mat3_mul(r1, r2);

    vecmath::row_mat3_mul(r1r2, r3)
}

/// Compute the rotation matrix for a set of 3-1-3 Euler angles (precession, nutation, spin)
pub fn euler313_to_dcm(theta: &Vector3<f64>) -> Matrix3<f64> {
    let r3 = euler3_rotation_degrees(theta[0]);
    let r2 = euler2_rotation_degrees(theta[1]);
    let r1 = euler3_rotation_degrees(theta[2]);

    let r1r2 = vecmath::row_mat3_mul(r1, r2);

    vecmath::row_mat3_mul(r1r2, r3)
}

/// Convert a quaternion to its corresponding rotation matrix
pub fn quat_to_dcm(q: &Quaternion<f64>) -> Matrix3<f64> {
    let q0 = q.0;
    let q1 = q.1[0];
    let q2 = q.1[1];
    let q3 = q.1[2];

    [
        [q0*q0 + q1*q1 - q2*q2 - q3*q3, 2.0*(q1*q2 - q0*q3), 2.0*(q1*q3 + q0*q2)],
        [2.0*(q1*q2 + q0*q3), q0*q0 - q1*q1 + q2*q2 - q3*q3, 2.0*(q2*q3 - q0*q1)],
        [2.0*(q1*q3 - q0*q2), 2.0*(q2*q3 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3],
    ]
}

/// Convert a rotation matrix to its corresponding quaternion
pub fn dcm_to_quat(dcm: &Matrix3<f64>) -> Quaternion<f64> {
    let trace = dcm[0][0] + dcm[1][1] + dcm[2][2];
    let q_squares = [
        (1.0 + trace)/4.0,
        (1.0 + 2.0*dcm[0][0] - trace)/4.0,
        (1.0 + 2.0*dcm[1][1] - trace)/4.0,
        (1.0 + 2.0*dcm[2][2] - trace)/4.0,
    ];
    let (max_index, max_q_square) = q_squares.iter().enumerate().max_by(|&(_, a), &(_, b)| a.partial_cmp(b).unwrap()).unwrap();
    let q0q1 = (dcm[1][2] - dcm[2][1]) / 4.0;
    let q0q2 = (dcm[2][0] - dcm[0][2]) / 4.0;
    let q0q3 = (dcm[0][1] - dcm[1][0]) / 4.0;
    let q1q2 = (dcm[0][1] + dcm[1][0]) / 4.0;
    let q1q3 = (dcm[0][2] + dcm[2][0]) / 4.0;
    let q2q3 = (dcm[1][2] + dcm[2][1]) / 4.0;

    let q: Quaternion<f64> = match max_index {
        0 => {
            let q0 = max_q_square.sqrt();
            (q0, [q0q1 / q0, q0q2 / q0, q0q3 / q0])
        },
        1 => {
            let q1 = max_q_square.sqrt();
            (q0q1 / q1, [q1, q1q2 / q1, q1q3 / q1])
        },
        2 => {
            let q2 = max_q_square.sqrt();
            (q0q2 / q2, [q1q2 / q2, q2, q2q3 / q2])
        },
        3 => {
            let q3 = max_q_square.sqrt();
            (q0q3 / q3, [q1q3 / q3, q2q3 / q3, q3])
        },
        _ => panic!("Invalid index"),
    };
    math::quat_normalize(q)
}

/// Differential Equation Matrix B for a quaternion
/// Such that qdot = 1/2 * B * [0, omega]'
/// where omega is the angular velocity,
/// and B is orthogonal :)
pub fn quat_differential_equation_matrix(q: Quaternion<f64>) -> Matrix4<f64> {
    let q0 = q.0;
    let q1 = q.1[0];
    let q2 = q.1[1];
    let q3 = q.1[2];

    [
        [q0, -q1, -q2, -q3],
        [q1, q0, -q3, q2],
        [q2, q3, q0, -q1],
        [q3, -q2, q1, q0],
    ]
}

/// Differential equation for quaternions, solve for qdot (time derivative of the quaternion) given omega (angular velocity)
pub fn omega_to_qdot(q: Quaternion<f64>, omega: Vector3<f64>) -> Quaternion<f64> {
    let b = quat_differential_equation_matrix(q);
    let v = [0.0, omega[0]*0.5, omega[1]*0.5, omega[2]*0.5];
    let qdot_vec = vecmath::row_mat4_transform(b, v);
    (qdot_vec[0], [qdot_vec[1], qdot_vec[2], qdot_vec[3]])
}

/// Differential equation for quaternions, solve for omega (angular velocity) given qdot (time derivative of the quaternion)
pub fn qdot_to_omega(q: Quaternion<f64>, q_dot: Quaternion<f64>) -> Vector3<f64> {
    let b = quat_differential_equation_matrix(q);
    let b_inv = vecmath::mat4_transposed(b); // B matrix is orthogonal!
    let q_dot_v = [q_dot.0*2.0, q_dot.1[0]*2.0, q_dot.1[1]*2.0, q_dot.1[2]*2.0];
    let omega_vec = vecmath::row_mat4_transform(b_inv, q_dot_v);
    [omega_vec[1], omega_vec[2], omega_vec[3]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dcm_to_quat() {
        // First with identity matrix:
        let dcm = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let q = dcm_to_quat(&dcm);
        assert!(math::quat_eq(&q, &(1.0, [0.0, 0.0, 0.0])));

        // Now with a rotation around y:
        let dcm = [
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let q = dcm_to_quat(&dcm);
        assert!(math::quat_eq(&q, &(0.0, [0.0, 1.0, 0.0])));
    }

    #[test]
    fn test_quat_to_dcm() {
        // First with identity quaternion:
        let q = (1.0, [0.0, 0.0, 0.0]);
        let dcm = quat_to_dcm(&q);
        assert!(math::mat3_eq(&dcm, &[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]));

        // Now with a rotation around y:
        let q = (0.0, [0.0, 1.0, 0.0]);
        let dcm = quat_to_dcm(&q);
        assert!(math::mat3_eq(&dcm, &[
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
        ]));
    }


}
