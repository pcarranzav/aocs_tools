//! Control module (TODO)

use vecmath::Vector3;

use crate::math::Vector6;

/// Convenience control or disturbance function that returns zero torque for any time and state vector
pub fn zero_torque(_: f64, _: &Vector6) -> Vector3<f64> {
    [0.0, 0.0, 0.0]
}
