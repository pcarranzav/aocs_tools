//! Plotting module

use crate::math::{Vector6, Vector7};
use gnuplot::*;

/// Plot the result of an MRP simulation, showing the MRP attitude and the angular velocity over time
/// This is a quick and dirty implementation that will open a new figure for each axis
pub fn plot_mrp_simulation(x: &Vec<Vector6>, t: &Vec<f64>) {
    let mut fg = Figure::new();
    let mrp0 = x.iter().map(|x| x[0]).collect::<Vec<_>>();
    let mrp1 = x.iter().map(|x| x[1]).collect::<Vec<_>>();
    let mrp2 = x.iter().map(|x| x[2]).collect::<Vec<_>>();
    let omega0 = x.iter().map(|x| x[3]).collect::<Vec<_>>();
    let omega1 = x.iter().map(|x| x[4]).collect::<Vec<_>>();
    let omega2 = x.iter().map(|x| x[5]).collect::<Vec<_>>();
    fg.axes2d()
        .set_title("Attitude (MRP X-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("MRP X-axis", &[])
        .lines(t, &mrp0, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude (MRP Y-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("MRP Y-axis", &[])
        .lines(t, &mrp1, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude (MRP Z-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("MRP Z-axis", &[])
        .lines(t, &mrp2, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate X", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate X-axis", &[])
        .lines(t, &omega0, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate Y", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate Y-axis", &[])
        .lines(t, &omega1, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate Z", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate Z-axis", &[])
        .lines(t, &omega2, &[]);
    fg.show_and_keep_running().unwrap();
}

/// Plot the result of an quaternion attitude simulation, showing the attitude and the angular velocity over time
/// This is a quick and dirty implementation that will open a new figure for each axis
pub fn plot_quat_simulation(x: &Vec<Vector7>, t: &Vec<f64>) {
    let mut fg = Figure::new();
    let q0 = x.iter().map(|x| x[0]).collect::<Vec<_>>();
    let q1 = x.iter().map(|x| x[1]).collect::<Vec<_>>();
    let q2 = x.iter().map(|x| x[2]).collect::<Vec<_>>();
    let q3 = x.iter().map(|x| x[3]).collect::<Vec<_>>();
    let omega0 = x.iter().map(|x| x[4]).collect::<Vec<_>>();
    let omega1 = x.iter().map(|x| x[5]).collect::<Vec<_>>();
    let omega2 = x.iter().map(|x| x[6]).collect::<Vec<_>>();
    fg.axes2d()
        .set_title("Attitude (Q Scalar)", &[])
        .set_x_label("t", &[])
        .set_y_label("Q Scalar", &[])
        .lines(t, &q0, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude (Q X-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("Q X-axis", &[])
        .lines(t, &q1, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude (Q Y-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("Q Y-axis", &[])
        .lines(t, &q2, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude (Q Z-axis)", &[])
        .set_x_label("t", &[])
        .set_y_label("Q Z-axis", &[])
        .lines(t, &q3, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate X", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate X-axis", &[])
        .lines(t, &omega0, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate Y", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate Y-axis", &[])
        .lines(t, &omega1, &[]);
    fg.show_and_keep_running().unwrap();
    let mut fg = Figure::new();
    fg.axes2d()
        .set_title("Attitude rate Z", &[])
        .set_x_label("t", &[])
        .set_y_label("Angular rate Z-axis", &[])
        .lines(t, &omega2, &[]);
    fg.show_and_keep_running().unwrap();
}
