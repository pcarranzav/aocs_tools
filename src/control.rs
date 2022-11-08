//! Control module (WIP)

use vecmath::Vector3;

use crate::math::{Vector6, Vector7};

/// Convenience control or disturbance function that returns zero torque for any time and MRP state vector
pub fn mrp_zero_torque(_: f64, _: &Vector6) -> Vector3<f64> {
    [0.0, 0.0, 0.0]
}

/// Convenience control or disturbance function that returns zero torque for any time and quaternion state vector
pub fn quat_zero_torque(_: f64, _: &Vector7) -> Vector3<f64> {
    [0.0, 0.0, 0.0]
}

/// Trait for controllers that work with MRP attitude representation
pub trait MRPController {
    /// For a time and state vector, return a control torque
    fn control_torque(&self, time: f64, state: &Vector6) -> Vector3<f64>;
}


/// Trait for controllers that work with quaternion attitude representation
pub trait QuatController {
    /// For a time and state vector, return a control torque
    fn control_torque(&self, time: f64, state: &Vector7) -> Vector3<f64>;
}

/// Controller for detumbling a rigid spacecraft using the MRP attitude representation
pub struct MRPDetumbleController {
    gain: f64,
}

impl MRPDetumbleController {
    /// Create a new detumble controller with a given gain
    pub fn new(gain: f64) -> Self {
        Self { gain }
    }
}

impl MRPController for MRPDetumbleController {
    fn control_torque(&self, _: f64, state: &Vector6) -> Vector3<f64> {
        [
            -self.gain * state[3],
            -self.gain * state[4],
            -self.gain * state[5],
        ]
    }
}

/// Controller for detumbling a rigid spacecraft using the quaternion attitude representation
pub struct QuatDetumbleController {
    gain: f64,
}

impl QuatDetumbleController {
    /// Create a new detumble controller with a given gain
    pub fn new(gain: f64) -> Self {
        Self { gain }
    }
}

impl QuatController for QuatDetumbleController {
    fn control_torque(&self, _: f64, state: &Vector7) -> Vector3<f64> {
        [
            -self.gain * state[4],
            -self.gain * state[5],
            -self.gain * state[6],
        ]
    }
}
