//! Simulation module (TODO)

use crate::math::VectorN;

type StateDerivativeFn<T> = fn(f64, &T) -> T;
type PostprocessingFn<T> = fn(&T) -> T;

/// Generic RK4 simulator.
/// It takes in the state derivative function, the initial conditions, the step size, the final time, and an optional postprocessing function
/// The postprocessing function sanitizes the state as in xn = p(xn); e.g. constraining its module or normalizing it.
/// The function returns the state and time vectors.
/// Note the state can be any size of vector, as it depends on the attitude representation or whatever we're trying to simulate,
/// the vector type must implement math::VectorN (addition and multiply by scalar).
pub fn rk4_simulator<T: VectorN>(f: StateDerivativeFn<T>, x0: T, step: f64, t_f: f64, p: Option<PostprocessingFn<T>>) -> (Vec<T>, Vec<f64>) {
    let mut t = 0.0;
    let mut x = vec![x0];
    let mut t_vec = vec![t];
    while t < t_f {
        let x_last = x.last().unwrap(); // Must exist because we're always including an element
        let k1 = f(t, &x_last).vec_scale(step);
        let t_mid = t + step / 2.0;
        let k2 = f(t_mid, &x_last.vec_add(&k1.vec_scale(0.5))).vec_scale(step);
        let k3 = f(t_mid, &x_last.vec_add(&k2.vec_scale(0.5))).vec_scale(step);
        t += step;
        let k4 = f(t, &x_last.vec_add(&k3)).vec_scale(step);

        let k2 = k2.vec_scale(2.0);
        let k3 = k3.vec_scale(2.0);

        let k_total = k1.vec_add(&k2).vec_add(&k3).vec_add(&k4).vec_scale(1.0 / 6.0);
        let mut xn = x_last.vec_add(&k_total);
        if let Some(p) = p {
            xn = p(&xn);
        }
        x.push(xn);
        t_vec.push(t);
    }
    (x, t_vec)
}
