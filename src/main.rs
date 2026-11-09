pub mod plot;

use crate::plot::{plot_convergence, plot_geometry, PlotSeries};
use complex_bessel::{besselj, hankel1};
use faer::Mat;
use faer::Side::Lower;
use num_complex::Complex64;
use std::error::Error;
use std::f64::consts::PI;

struct Panel {
    a: f64,
    b: f64,
    points: Vec<Complex64>,
    weights: Vec<f64>,
    densities: Vec<Complex64>,
    length: f64,
}

struct Boundary {
    curve: Curve,
    panels: Vec<Panel>,
}

#[derive(Clone, Copy)]
enum Curve {
    Circle,
    Ellipse { alpha: f64 },
    Starfish,
}

impl Curve {
    fn x(&self, t: f64) -> Complex64 {
        let theta = 2.0 * PI * t;
        match self {
            Curve::Circle => Complex64::new(theta.cos(), theta.sin()),
            Curve::Ellipse { alpha } => Complex64::new(theta.cos(), theta.sin() / alpha),
            Curve::Starfish => {
                let r = 1.0 + 0.25 * (5.0 * theta).sin();
                Complex64::from_polar(r, theta)
            }
        }
    }

    fn dx(&self, t: f64) -> Complex64 {
        let theta = 2.0 * PI * t;
        match self {
            Curve::Circle => Complex64::new(-2.0 * PI * theta.sin(), 2.0 * PI * theta.cos()),
            Curve::Ellipse { alpha } => Complex64::new(-2.0 * PI * theta.sin(), 2.0 * PI * theta.cos() / alpha),
            Curve::Starfish => {
                let r = 1.0 + 0.25 * (5.0 * theta).sin();
                let dr = 2.5 * PI * (5.0 * theta).cos();
                Complex64::from_polar(1.0, theta) * Complex64::new(dr, 2.0 * PI * r)
            }
        }
    }

    fn normal(&self, t: f64) -> Complex64 {
        let tangent = self.dx(t);
        Complex64::new(tangent.im, -tangent.re) / tangent.norm()
    }
}

fn density(t: f64, m: i32) -> Complex64 { // fourier for easy analysis
    Complex64::from_polar(1.0, 2.0 * PI * m as f64 * t)
}

fn gauss_legendre_nodes_weights(n: usize) -> Result<(Vec<f64>, Vec<f64>), Box<dyn Error>> {
    let J = Mat::from_fn(n, n, |i, j| {
        if i == j + 1 {
            (i as f64) / (4.0 * (i as f64).powf(2.0) - 1.0).sqrt()
        } else {
            0.0
        }
    });
    let eigen = J.self_adjoint_eigen(Lower).unwrap(); // evderror hello??
    Ok((eigen.S().column_vector().iter().map(|x| *x).collect(),  (0..n).map(|j| 2.0 * eigen.U()[(0, j)].powi(2)).collect()))
}

fn legendre(n: usize, x: f64) -> f64 {
    if n == 0 {
        1.0
    } else if n == 1 {
        x
    } else {
        let mut a = 1.0;
        let mut b = x;
        for i in 1..n {
            let c = ((2 * i + 1) as f64 * x * b - i as f64 * a) / (i + 1) as f64;
            a = b;
            b = c;
        }
        b
    }
}


fn panel_length(curve: Curve, a: f64, b: f64, nodes: &[f64], weights: &[f64]) -> f64 {
    nodes.iter()
        .zip(weights.iter())
        .map(|(&s, &w)| {
            let t = (a + b) / 2.0 + (b - a) / 2.0 * s;
            (b - a) / 2.0 * w * curve.dx(t).norm()
        }).sum()
}

impl Boundary {
    fn new(curve: Curve, M: usize, q: usize, m: i32) -> Self {
        let (nodes, weights) = gauss_legendre_nodes_weights(64).unwrap();
        let panels: Vec<(f64, f64)> = (0..M).map(|i| ((i as f64 ) / M as f64, (i + 1) as f64 / M as f64)).collect();
        let mut panels: Vec<(f64, f64, f64)> = panels.iter().map(|(a, b)| (*a, *b, panel_length(curve, *a, *b, &nodes, &weights))).collect();

        loop { // naive bullshit but im tired and 50 is a small number
            let mut new_panels = Vec::new();
            for (i, (a, b, length)) in panels.iter().enumerate() {
                let midpoint = (a + b) / 2.0;
                let prev_length = panels[if i != 0 { i - 1 } else { panels.len() - 1 }].2;
                let next_length = panels[if i != (panels.len() - 1) { i + 1 } else { 0 }].2;

                if *length > 2.0 * prev_length || *length > 2.0 * next_length {
                    new_panels.push((*a, midpoint, panel_length(curve, *a, midpoint, &nodes, &weights)));
                    new_panels.push((midpoint, *b, panel_length(curve, midpoint, *b, &nodes, &weights)));
                    continue;
                }

                let mut tail_energy = 0.0;

                for i in 16..64 {
                    let coefficient = (2 * i + 1) as f64 / 2.0 * nodes.iter().zip(weights.iter()).map(|(&s, &w)| {
                        let t = midpoint + (b - a) / 2.0 * s;
                        w * curve.dx(t).norm() * legendre(i, s)
                    }).sum::<f64>();
                    tail_energy += 2.0 / (2 * i + 1) as f64 * coefficient.powf(2.0);
                }
                tail_energy = tail_energy.sqrt();

                if tail_energy > 1e-11 {
                    new_panels.push((*a, midpoint, panel_length(curve, *a, midpoint, &nodes, &weights)));
                    new_panels.push((midpoint, *b, panel_length(curve, midpoint, *b, &nodes, &weights)));
                } else {
                    new_panels.push((*a, *b, *length));
                }
            }
            if new_panels.len() == panels.len() {
                break;
            }
            panels = new_panels;
        }

        let (nodes, weights) = gauss_legendre_nodes_weights(q).unwrap();
        let mut p = Vec::new();

        for (a, b, length) in panels {
            let (points, weights, densities): (Vec<Complex64>, Vec<f64>, Vec<Complex64>) = nodes.iter().zip(weights.iter()).enumerate().map(|(j, (&node, &weight))| {
                let t = 0.5 * (b - a) * node + 0.5 * (b + a);
                let x = curve.x(t);
                (x, 0.5 * (b - a) * weight * curve.dx(t).norm(), density(t, m))
            }).collect();
            let length = weights.iter().sum();
            p.push(Panel {
                a,
                b,
                points,
                weights,
                densities,
                length,
            });
        }

        Self { curve, panels: p }
    }
}

fn qbx(boundary: &Boundary, t_target: f64, k: f64, p: i32, radius_factor: f64) -> Result<Complex64, Box<dyn Error>> {
    let panel_index = boundary.panels
        .iter()
        .position(|panel| {
            t_target >= panel.a && t_target < panel.b
    }).unwrap();
    let h = boundary.panels[panel_index].length;
    let radius =  h * radius_factor;
    let target = boundary.curve.x(t_target);
    let center = target - boundary.curve.normal(t_target) * radius;

    let mut alpha = Vec::new();
    for l in -p..=p {
        let mut sum = Complex64::new(0.0, 0.0);
        for panel in boundary.panels.iter() {
            for i in 0..panel.points.len() {
                let d = panel.points[i] - center;
                sum += hankel1(l as f64, Complex64::new(k * d.norm(), 0.0))? * Complex64::from_polar(1.0, l as f64 * d.arg()) * panel.densities[i] * panel.weights[i];
            }
        }
        alpha.push(Complex64::new(0.0, 0.25) * sum);
    }

    let mut u = Complex64::ZERO;
    for l in -p..=p {
        let d = target - center;
        u += alpha[(l + p) as usize] * besselj(l as f64, Complex64::new(k * d.norm(), 0.0))? * Complex64::from_polar(1.0, -l as f64 * d.arg());
    }
    Ok(u)
}

fn exact_circle_solution(k: f64, m: i32, target: Complex64) -> Result<Complex64, Box<dyn Error>> {
    let k_boundary = Complex64::new(k, 0.0);
    let k_target = Complex64::new(k * target.norm(), 0.0);

    let radial = if target.norm() <= 1.0 {
        hankel1(m as f64, k_boundary)? * besselj(m as f64, k_target)?
    } else {
        besselj(m as f64, k_boundary)? * hankel1(m as f64, k_target)?
    };

    Ok(Complex64::new(0.0, PI / 2.0) * radial * Complex64::from_polar(1.0, m as f64 * target.arg()))
}

fn evaluate_params(curve: Curve, k: f64, p: i32, M: usize, q: usize, m: i32, radius_factor: f64, t_target: f64) -> Result<(Complex64, Complex64), Box<dyn Error>> {
    let boundary = Boundary::new(curve, M, q, m);
    let target = curve.x(t_target);

    let u = qbx(&boundary, t_target, k, p, radius_factor)?;
    let exact = exact_circle_solution(k, m, target)?;

    Ok((u, exact))
}


fn max_reference_error(boundary: &Boundary, reference: &Boundary, targets: &[f64], k: f64, p: i32, reference_p: i32) -> Result<f64, Box<dyn Error>> {
    let mut max_error: f64 = 0.0;

    for &t in targets {
        let computed = qbx(boundary, t, k, p, 0.5)?;
        let expected = qbx(reference, t, k, reference_p, 0.5)?;
        max_error = max_error.max((computed - expected).norm());
    }

    Ok(max_error)
}


fn noncircular_convergence() -> Result<(), Box<dyn Error>> {
    let k = 1.0;
    let m = 2;
    let p = 8;
    let q = 32;

    let targets: Vec<f64> = (0..12)
        .map(|j| (j as f64 + 0.37) / 12.0)
        .collect();

    let geometries = [
        ("Ellipse alpha=3", Curve::Ellipse { alpha: 3.0 }),
        ("Ellipse alpha=6", Curve::Ellipse { alpha: 6.0 }),
        ("Starfish", Curve::Starfish),
    ];

    let mut series = Vec::new();

    for (name, curve) in geometries {
        let reference = Boundary::new(curve, 200, 96, m);

        let mut points = Vec::new();

        for initial_m in [10, 20, 40, 80, 160] {
            let boundary = Boundary::new(
                curve,
                initial_m,
                q,
                m,
            );

            let error = max_reference_error(
                &boundary,
                &reference,
                &targets,
                k,
                p,
                16,
            )?;
            points.push((initial_m as f64, error));
        }

        series.push(PlotSeries {
            label: name.to_string(),
            points,
        });
    }

    plot_convergence(
        "plots/noncircular_convergence.svg",
        "QBX Convergence on Noncircular Geometries",
        "Initial panel count M",
        &series,
    )?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    std::fs::create_dir_all("plots")?;
    let k = 1.0; // wave number
    let p = 8; // expansion order
    let M = 50; // panels
    let m = 2; // mode
    let radius_factor = 0.5;
    let t_target = 0.25;

    let mut series = Vec::new();

    for q in [4, 8, 16, 32, 64] {
        let mut points = Vec::new();

        for p in [2, 4, 6, 8, 10, 12, 16, 20] {
            let (u, exact) = evaluate_params(Curve::Circle, k, p, M, q, m, radius_factor, t_target)?;
            let error = (u - exact).norm();
            points.push((p as f64, error));
        }

        series.push(PlotSeries {
            label: format!("q = {q}"),
            points,
        });
    }

    plot_convergence(
        "plots/error_vs_order.svg",
        "QBX Error vs Expansion Order",
        "Expansion order p",
        &series,
    )?;

    let mut series = Vec::new();

    let radius_factors = [
        0.25, 0.375, 0.5, 0.625,
        0.75, 0.875, 1.0,
    ];

    for q in [8, 16, 32, 64] {
        let mut points = Vec::new();

        for radius_factor in radius_factors {
            let (u, exact) = evaluate_params(Curve::Circle, k, p, M, q, m, radius_factor, t_target)?;
            let error = (u - exact).norm();
            points.push((radius_factor, error));
        }

        series.push(PlotSeries {
            label: format!("q = {q}"),
            points,
        });
    }

    plot_convergence(
        "plots/error_vs_radius.svg",
        "QBX Error vs Center Distance",
        "Radius ratio r/h",
        &series,
    )?;

    let boundary = Boundary::new(Curve::Starfish, 50, 16, 2);

    plot_geometry(
        "plots/starfish_geometry.svg",
        &boundary,
        0.5,
    )?;

    noncircular_convergence()?;

    Ok(())
}
