//! Simulation module (TODO)

use crate::math;

type StateDerivativeFn = fn(f64, &Vec<f64>) -> Vec<f64>;
type PostprocessingFn = fn(&Vec<f64>) -> Vec<f64>;

/// Generic RK4 simulator.
/// It takes in the state derivative function, the initial conditions, the step size, the final time, and an optional postprocessing function
/// The postprocessing function sanitizes the state as in xn = p(xn); e.g. constraining its module or normalizing it.
/// The function returns the state and time vectors.
/// Note the state can be any size of vector, as it depends on the attitude representation or whatever we're trying to simulate.
pub fn rk4_simulator(f: StateDerivativeFn, x0: Vec<f64>, step: f64, t_f: f64, p: Option<PostprocessingFn>) -> (Vec<Vec<f64>>, Vec<f64>) {
    let mut t = 0.0;
    let mut x = vec![x0];
    let mut t_vec = vec![t];
    while t < t_f {
        let x_last = x.last().unwrap(); // Must exist because we're always including an element
        let k1 = math::vecn_scale(&f(t, &x_last), step);
        let t_mid = t + step / 2.0;
        let k2 = math::vecn_scale(&f(t_mid, &math::vecn_add(&x_last, &math::vecn_scale(&k1, 0.5))), step);
        let k3 = math::vecn_scale(&f(t_mid, &math::vecn_add(&x_last, &math::vecn_scale(&k2, 0.5))), step);
        t += step;
        let k4 = math::vecn_scale(&f(t, &math::vecn_add(&x_last, &k3)), step);

        let k2 = math::vecn_scale(&k2, 2.0);
        let k3 = math::vecn_scale(&k3, 2.0);

        let k_total = math::vecn_add(&math::vecn_add(&math::vecn_add(&k1, &k2), &k3), &k4);
        let k_total = math::vecn_scale(&k_total, 1.0/6.0);
        let mut xn = math::vecn_add(&x_last, &k_total);
        if let Some(p) = p {
            xn = p(&xn);
        }
        x.push(xn);
        t_vec.push(t);
    }
    (x, t_vec)
}
