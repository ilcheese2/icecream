use plotters::prelude::*; // this file was entirely created by ai
use std::error::Error;
use crate::Boundary;

pub struct PlotSeries {
    pub label: String,
    pub points: Vec<(f64, f64)>,
}

pub fn plot_convergence(
    filename: &str,
    title: &str,
    x_label: &str,
    series: &[PlotSeries],
) -> Result<(), Box<dyn Error>> {
    // Convert errors to log10(error).
    // Filter out invalid/zero values.
    let transformed: Vec<Vec<(f64, f64)>> = series
        .iter()
        .map(|s| {
            s.points
                .iter()
                .filter(|(x, e)| {
                    x.is_finite() && e.is_finite() && *e > 0.0
                })
                .map(|(x, e)| (*x, e.log10()))
                .collect()
        })
        .collect();

    let all: Vec<(f64, f64)> = transformed
        .iter()
        .flat_map(|s| s.iter().copied())
        .collect();

    if all.is_empty() {
        return Err("No valid data to plot".into());
    }

    let x_min = all.iter().map(|p| p.0)
        .fold(f64::INFINITY, f64::min);
    let x_max = all.iter().map(|p| p.0)
        .fold(f64::NEG_INFINITY, f64::max);

    let y_min = all.iter().map(|p| p.1)
        .fold(f64::INFINITY, f64::min);
    let y_max = all.iter().map(|p| p.1)
        .fold(f64::NEG_INFINITY, f64::max);

    let x_pad = ((x_max - x_min).abs() * 0.05).max(0.01);
    let y_pad = ((y_max - y_min).abs() * 0.10).max(0.5);

    let root = SVGBackend::new(filename, (1100, 750))
        .into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 32))
        .margin(25)
        .x_label_area_size(65)
        .y_label_area_size(85)
        .build_cartesian_2d(
            (x_min - x_pad)..(x_max + x_pad),
            (y_min - y_pad)..(y_max + y_pad),
        )?;

    chart.configure_mesh()
        .x_desc(x_label)
        .y_desc("log10(absolute error)")
        .axis_desc_style(("sans-serif", 22))
        .label_style(("sans-serif", 16))
        .draw()?;

    for (i, data) in transformed.iter().enumerate() {
        if data.is_empty() {
            continue;
        }

        let color = Palette99::pick(i);

        chart
            .draw_series(LineSeries::new(
                data.iter().copied(),
                color.stroke_width(3),
            ))?
            .label(series[i].label.clone())
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 25, y)],
                    Palette99::pick(i).stroke_width(3),
                )
            });

        chart.draw_series(
            data.iter().map(|&(x, y)| {
                Circle::new(
                    (x, y),
                    4,
                    Palette99::pick(i).filled(),
                )
            }),
        )?;
    }

    chart.configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .label_font(("sans-serif", 18))
        .draw()?;

    root.present()?;
    Ok(())
}

pub(crate) fn plot_geometry(
    filename: &str,
    boundary: &Boundary,
    radius_factor: f64,
) -> Result<(), Box<dyn Error>> {
    let root = SVGBackend::new(filename, (900, 900))
        .into_drawing_area();

    root.fill(&WHITE)?;

    let samples: Vec<(f64, f64)> = (0..=1000)
        .map(|i| {
            let t = i as f64 / 1000.0;
            let z = boundary.curve.x(t);
            (z.re, z.im)
        })
        .collect();

    let mut chart = ChartBuilder::on(&root)
        .caption("QBX Expansion Geometry", ("sans-serif", 28))
        .margin(20)
        .build_cartesian_2d(-1.8..1.8, -1.8..1.8)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
        samples,
        BLUE.stroke_width(3),
    ))?;

    let n = boundary.panels.len();

    // Plot one center/disk per panel at its midpoint.
    for i in 0..n {
        let panel = &boundary.panels[i];
        let t = 0.5 * (panel.a + panel.b);

        let target = boundary.curve.x(t);
        let normal = boundary.curve.normal(t);

        let radius = boundary.panels[i].length * radius_factor;
        let center = target - radius * normal;

        // Expansion circle in physical coordinates
        let disk = (0..=80).map(|j| {
            let angle = 2.0 * std::f64::consts::PI
                * j as f64 / 80.0;

            (
                center.re + radius * angle.cos(),
                center.im + radius * angle.sin(),
            )
        });

        chart.draw_series(LineSeries::new(
            disk,
            GREEN.mix(0.3).stroke_width(1),
        ))?;

        chart.draw_series(std::iter::once(
            Circle::new(
                (center.re, center.im),
                2,
                RED.filled(),
            )
        ))?;
    }

    root.present()?;
    Ok(())
}