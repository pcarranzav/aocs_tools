# AOCS tools (WIP)

This crate contains tools for Attitude and Orbit Control Systems.

I'm mostly writing this to practice as I learn Rust, so do not treat it as production-ready code in any way.

It is early work in progress, but the plan is to include utilities for:

- Kinematics: computing and describing rotations using various attitude representations (Euler angles, rotation matrices, quaternions, Modified Rodrigues parameters).
- Estimation (TODO): Producing attitude values based on sensor readings. I'll include Triad, QUEST and OLAE probably. I may or may not add an Extended Kalman Filter, depending on whether I have time for it...
- Control: applying feedback laws to modify the attitude and/or rate to match a desired target. For now I only included basic detumbling, but I'd like to add more controllers with various targeting modes.
- Simulation: using numerical methods to analyze specific instances of a body rotating based on its kinetics and including a control law. For now we include a basic Runge-Kutta 4 simulator, with implementations for quaternion and MRP based state representations. The simulation accepts a function to define the disturbance torques, but for now we always set this to zero; it'd be great to add some models for drag, magnetic dipole, etc.

This is mostly based on the courses provided in the [Spacecraft Dynamics and Control Specialization in Coursera](https://www.coursera.org/specializations/spacecraft-dynamics-control?), by Professor Schaub from CU-Boulder. We'll use Professor Schaub's notation.

In particular, note that for a rotation matrix between a body frame B and an inertial frame N, we'll call this matrix BN to represent:

```
// vector_in_B_frame = BN * vector_in_N_frame
```

i.e. the matrix that, multiplied by a vector in the inertial frame, gives us the vector in the body frame.

We call rotation matrices "DCM" i.e. the Direction Cosine Matrix.

For quaternions, we use the [quaternion](https://docs.rs/quaternion/latest/quaternion/) crate, which uses the following definition of the quaternion type:

```rust
pub type Quaternion<T> = (T, [T; 3]);
```

This means the scalar part of the quaternion will be the first item in the tuple, and the second item contains the vector part.
