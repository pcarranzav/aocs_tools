//! Kinematics module (WIP)
//! Provides tools to describe rotations/attitudes and to
//! convert between different representations of rotations.

use quaternion::Quaternion;
use vecmath::{Vector3, Matrix3};

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
pub fn quaternion_to_dcm(q: &Quaternion<f64>) -> Matrix3<f64> {
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
pub fn dcm_to_quaternion(dcm: &Matrix3<f64>) -> Quaternion<f64> {
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
    math::quaternion_normalize(q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dcm_to_quaternion() {
        // First with identity matrix:
        let dcm = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let q = dcm_to_quaternion(&dcm);
        assert!(math::quaternion_eq(&q, &(1.0, [0.0, 0.0, 0.0])));

        // Now with a rotation around y:
        let dcm = [
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let q = dcm_to_quaternion(&dcm);
        assert!(math::quaternion_eq(&q, &(0.0, [0.0, 1.0, 0.0])));
    }

    #[test]
    fn test_quaternion_to_dcm() {
        // First with identity quaternion:
        let q = (1.0, [0.0, 0.0, 0.0]);
        let dcm = quaternion_to_dcm(&q);
        assert!(math::mat3_eq(&dcm, &[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]));

        // Now with a rotation around y:
        let q = (0.0, [0.0, 1.0, 0.0]);
        let dcm = quaternion_to_dcm(&q);
        assert!(math::mat3_eq(&dcm, &[
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
        ]));
    }


}
