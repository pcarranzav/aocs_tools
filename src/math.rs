//! Math utilities

use std::f64::EPSILON;

use quaternion::Quaternion;
use vecmath::Matrix3;
use vecmath::Vector3;

/// Good old Pi
pub const PI: f64 = 3.14159265358979323846264338327950288_f64;


/// 6D float vector, useful for simulations using MRPs
pub type Vector6 = [f64; 6];
/// 7D float vector, useful for simulations using quaternions
pub type Vector7 = [f64; 7];
/// Trait for arrays behaving as arbitrary sized vectors
pub trait VectorN {
    /// Add the vector to another vector
    fn vec_add(&self, other: &Self) -> Self;
    /// Scale the vector by a scalar
    fn vec_scale(&self, scalar: f64) -> Self;
}

impl VectorN for [f64; 6] {
    fn vec_add(&self, other: &Self) -> Self {
        [
            self[0] + other[0],
            self[1] + other[1],
            self[2] + other[2],
            self[3] + other[3],
            self[4] + other[4],
            self[5] + other[5],
        ]
    }

    fn vec_scale(&self, scalar: f64) -> Self {
        [
            self[0] * scalar,
            self[1] * scalar,
            self[2] * scalar,
            self[3] * scalar,
            self[4] * scalar,
            self[5] * scalar,
        ]
    }
}

impl VectorN for [f64; 7] {
    fn vec_add(&self, other: &Self) -> Self {
        [
            self[0] + other[0],
            self[1] + other[1],
            self[2] + other[2],
            self[3] + other[3],
            self[4] + other[4],
            self[5] + other[5],
            self[6] + other[6],
        ]
    }

    fn vec_scale(&self, scalar: f64) -> Self {
        [
            self[0] * scalar,
            self[1] * scalar,
            self[2] * scalar,
            self[3] * scalar,
            self[4] * scalar,
            self[5] * scalar,
            self[6] * scalar,
        ]
    }
}

/// Compute the tilde matrix for a vector (i.e. matrix representation of the cross product)
pub fn vec3_tilde(vec: Vector3<f64>) -> Matrix3<f64> {
    [
        [0.0, -vec[2], vec[1]],
        [vec[2], 0.0, -vec[0]],
        [-vec[1], vec[0], 0.0],
    ]
}

/// Compute the outer product between two 3D vectors
pub fn vec3_outer_product(a: Vector3<f64>, b: Vector3<f64>) -> Matrix3<f64> {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}

/// Scale a vector of arbitrary size by a scalar
pub fn vecn_scale(vec: &Vec<f64>, scalar: f64) -> Vec<f64> {
    vec.iter().map(|&x| x * scalar).collect()
}

/// Add two vectors of arbitrary size
pub fn vecn_add(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Multiply a matrix by a scalar
pub fn mat3_scale(mat: Matrix3<f64>, scale: f64) -> Matrix3<f64> {
    [
        [mat[0][0]*scale, mat[0][1]*scale, mat[0][2]*scale],
        [mat[1][0]*scale, mat[1][1]*scale, mat[1][2]*scale],
        [mat[2][0]*scale, mat[2][1]*scale, mat[2][2]*scale],
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

/// Compute whether two 3D vectors are equal within float precision
pub fn vec3_eq(a: &Vector3<f64>, b: &Vector3<f64>) -> bool {
    for i in 0..3 {
        if (a[i] - b[i]).abs() > EPSILON {
            return false;
        }
    }
    true
}

/// Compute whether two 3x3 matrices are equal within float precision
pub fn mat3_eq(a: &Matrix3<f64>, b: &Matrix3<f64>) -> bool {
    for i in 0..3 {
        if !vec3_eq(&a[i], &b[i]) {
            return false;
        }
    }
    true
}

/// Normalize a quaternion so that the scalar part is positive
pub fn quat_positive_scalar(q: Quaternion<f64>) -> Quaternion<f64> {
    if q.0 < 0.0 {
        (-q.0, [-q.1[0], -q.1[1], -q.1[2]])
    } else {
        q
    }
}

/// Normalize a quaternion so that it has unit length
pub fn quat_unit_length(q: Quaternion<f64>) -> Quaternion<f64> {
    let len = quaternion::len(q);
    (q.0 / len, [q.1[0] / len, q.1[1] / len, q.1[2] / len])
}

/// Normalize a quaternion so that it has unit length and the scalar part is positive
pub fn quat_normalize(q: Quaternion<f64>) -> Quaternion<f64> {
    quat_positive_scalar(quat_unit_length(q))
}

/// Computes whether two (already normalized) quaternions are equal within float precision
pub fn quat_eq(a: &Quaternion<f64>, b: &Quaternion<f64>) -> bool {
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
    fn test_quat_normalize() {
        // Note the norm is 8.0:
        let q = quat_normalize((-4.0, [2.0, -6.0, (8.0_f64).sqrt()]));
        assert!((quaternion::len(q) - 1.0).abs() < EPSILON);
        assert!(quat_eq(&q, &(0.5, [-0.25, 0.75, -1.0/(8.0_f64).sqrt()])));
    }
}
