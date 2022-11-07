//! Math utilities

use std::f64::EPSILON;

use quaternion::Quaternion;
use vecmath::Matrix3;
use vecmath::Vector3;

/// Good old Pi
pub const PI: f64 = 3.14159265358979323846264338327950288_f64;

/// Compute the tilde matrix for a vector
pub fn vec3_tilde(vec: Vector3<f64>) -> Matrix3<f64>
{
    [
        [0.0, -vec[2], vec[1]],
        [vec[2], 0.0, -vec[0]],
        [-vec[1], vec[0], 0.0],
    ]
}

/// Compute the cosine of a number in degrees
pub fn cosd(x: f64) -> f64 {
    x.to_radians().cos()
}

/// Compute the sine of a number in degrees
pub fn sind(x: f64) -> f64 {
    x.to_radians().sin()
}

/// Compute whether two 3x3 matrices are equal within float precision
pub fn mat3_eq(a: &Matrix3<f64>, b: &Matrix3<f64>) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if (a[i][j] - b[i][j]).abs() > EPSILON {
                return false;
            }
        }
    }
    true
}

/// Normalize a quaternion so that the scalar part is positive
pub fn quaternion_positive_scalar(q: Quaternion<f64>) -> Quaternion<f64> {
    if q.0 < 0.0 {
        (-q.0, [-q.1[0], -q.1[1], -q.1[2]])
    } else {
        q
    }
}

/// Normalize a quaternion so that it has unit length
pub fn quaternion_unit_length(q: Quaternion<f64>) -> Quaternion<f64> {
    let len = quaternion::len(q);
    (q.0 / len, [q.1[0] / len, q.1[1] / len, q.1[2] / len])
}

/// Normalize a quaternion so that it has unit length and the scalar part is positive
pub fn quaternion_normalize(q: Quaternion<f64>) -> Quaternion<f64> {
    quaternion_positive_scalar(quaternion_unit_length(q))
}

/// Computes whether two (already normalized) quaternions are equal within float precision
pub fn quaternion_eq(a: &Quaternion<f64>, b: &Quaternion<f64>) -> bool {
    for i in 0..3 {
        if (a.1[i] - b.1[i]).abs() > EPSILON {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tilde() {
        let result: Matrix3<f64> = vec3_tilde([1.0, 2.0, 3.0]);
        let expected:  Matrix3<f64> = [
            [0.0, -3.0, 2.0],
            [3.0, 0.0, -1.0],
            [-2.0, 1.0, 0.0],
        ];
        assert!(mat3_eq(&result, &expected));
    }

    #[test]
    fn test_quaternion_normalize() {
        // Note the norm is 8.0:
        let q = quaternion_normalize((-4.0, [2.0, -6.0, (8.0_f64).sqrt()]));
        assert!((quaternion::len(q) - 1.0).abs() < EPSILON);
        assert!(quaternion_eq(&q, &(0.5, [-0.25, 0.75, -1.0/(8.0_f64).sqrt()])));
    }
}
