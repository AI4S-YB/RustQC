//! SVG plot generation for ATAC-seq QC metrics.
//!
//! Three plots: fragment size distribution, TSSE profile, library complexity saturation curve.

use std::path::Path;

use anyhow::Result;
use plotters::prelude::*;
use plotters_svg::SVGBackend;

use crate::atac::lib_complexity::LibComplexityRow;

// ---------------------------------------------------------------------------
// Fragment size distribution plot
// ---------------------------------------------------------------------------

/// Render a fragment-size density line plot to an SVG file.
///
/// X-axis: fragment length 0–1010 bp.
/// Y-axis: density × 1000 (per-mille).
pub fn fragsize_svg(rows: &[(u32, u64, f64)], path: &Path, sample: &str) -> Result<()> {
    let root = SVGBackend::new(path, (640, 400)).into_drawing_area();
    root.fill(&WHITE)?;

    let (top, plot_area) = root.split_vertically(40u32);
    let cx = 640 / 2;
    top.draw(&Text::new(
        "Fragment Size Distribution",
        (cx, 4),
        ("sans-serif", 16u32)
            .into_font()
            .style(FontStyle::Bold)
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;
    top.draw(&Text::new(
        sample,
        (cx, 22),
        ("sans-serif", 12u32)
            .into_font()
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;

    // Only draw non-zero region to determine y_max.
    let data: Vec<(f64, f64)> = rows
        .iter()
        .map(|&(l, _, d)| (l as f64, d * 1000.0))
        .collect();

    let y_max = data
        .iter()
        .map(|(_, y)| *y)
        .fold(0.0f64, f64::max)
        .max(1.0);

    let mut chart = ChartBuilder::on(&plot_area)
        .margin_top(5u32)
        .margin_right(20u32)
        .margin_bottom(10u32)
        .margin_left(10u32)
        .x_label_area_size(40u32)
        .y_label_area_size(55u32)
        .build_cartesian_2d(0f64..1010f64, 0f64..y_max * 1.05)?;

    chart
        .configure_mesh()
        .bold_line_style(TRANSPARENT)
        .light_line_style(TRANSPARENT)
        .x_desc("Fragment length (bp)")
        .y_desc("Density ×1000")
        .axis_desc_style(("sans-serif", 12u32))
        .label_style(("sans-serif", 10u32))
        .draw()?;

    chart.draw_series(LineSeries::new(data, BLUE.stroke_width(1)))?;

    plot_area.present()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// TSSE profile plot
// ---------------------------------------------------------------------------

/// Render a TSSE window-profile line plot to an SVG file.
///
/// X-axis: window index 1–20.
/// Y-axis: normalised signal.
pub fn tsse_svg(values: &[f64], path: &Path, sample: &str) -> Result<()> {
    let root = SVGBackend::new(path, (480, 360)).into_drawing_area();
    root.fill(&WHITE)?;

    let (top, plot_area) = root.split_vertically(40u32);
    let cx = 480 / 2;
    top.draw(&Text::new(
        "TSSE Score Profile",
        (cx, 4),
        ("sans-serif", 16u32)
            .into_font()
            .style(FontStyle::Bold)
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;
    top.draw(&Text::new(
        sample,
        (cx, 22),
        ("sans-serif", 12u32)
            .into_font()
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;

    let n = values.len();
    let data: Vec<(f64, f64)> = values
        .iter()
        .enumerate()
        .map(|(i, &v)| ((i + 1) as f64, v))
        .collect();

    let y_max = values
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .fold(0.0f64, f64::max)
        .max(1.0);

    let x_max = (n as f64).max(20.0) + 0.5;

    let mut chart = ChartBuilder::on(&plot_area)
        .margin_top(5u32)
        .margin_right(20u32)
        .margin_bottom(10u32)
        .margin_left(10u32)
        .x_label_area_size(40u32)
        .y_label_area_size(50u32)
        .build_cartesian_2d(0.5f64..x_max, 0f64..y_max * 1.1)?;

    chart
        .configure_mesh()
        .bold_line_style(TRANSPARENT)
        .light_line_style(TRANSPARENT)
        .x_desc("Window index")
        .y_desc("Normalised signal")
        .axis_desc_style(("sans-serif", 12u32))
        .label_style(("sans-serif", 10u32))
        .draw()?;

    chart.draw_series(LineSeries::new(data, RED.stroke_width(2)))?;

    // Mark the peak.
    if let Some((peak_i, &peak_v)) = values
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite())
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    {
        chart.draw_series(std::iter::once(Circle::new(
            ((peak_i + 1) as f64, peak_v),
            4u32,
            RED.filled(),
        )))?;
    }

    plot_area.present()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Library complexity saturation curve
// ---------------------------------------------------------------------------

/// Render a library complexity saturation curve to an SVG file.
///
/// X-axis: putative reads (millions).
/// Y-axis: distinct fragments (millions).
pub fn lib_complexity_svg(rows: &[LibComplexityRow], path: &Path, sample: &str) -> Result<()> {
    let root = SVGBackend::new(path, (480, 360)).into_drawing_area();
    root.fill(&WHITE)?;

    let (top, plot_area) = root.split_vertically(40u32);
    let cx = 480 / 2;
    top.draw(&Text::new(
        "Library Complexity",
        (cx, 4),
        ("sans-serif", 16u32)
            .into_font()
            .style(FontStyle::Bold)
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;
    top.draw(&Text::new(
        sample,
        (cx, 22),
        ("sans-serif", 12u32)
            .into_font()
            .color(&BLACK)
            .pos(plotters::style::text_anchor::Pos::new(
                plotters::style::text_anchor::HPos::Center,
                plotters::style::text_anchor::VPos::Top,
            )),
    ))?;

    // Convert to millions; skip NaN rows.
    let data: Vec<(f64, f64)> = rows
        .iter()
        .filter(|r| r.distinct_fragments.is_finite())
        .map(|r| (r.putative_reads / 1e6, r.distinct_fragments / 1e6))
        .collect();

    if data.is_empty() {
        // Write an empty SVG with a note.
        top.draw(&Text::new(
            "(no data)",
            (cx, 80),
            ("sans-serif", 14u32).into_font().color(&BLACK).pos(
                plotters::style::text_anchor::Pos::new(
                    plotters::style::text_anchor::HPos::Center,
                    plotters::style::text_anchor::VPos::Top,
                ),
            ),
        ))?;
        plot_area.present()?;
        return Ok(());
    }

    let x_max = data.iter().map(|(x, _)| *x).fold(0.0f64, f64::max).max(1.0) * 1.05;
    let y_max = data.iter().map(|(_, y)| *y).fold(0.0f64, f64::max).max(1.0) * 1.05;

    let mut chart = ChartBuilder::on(&plot_area)
        .margin_top(5u32)
        .margin_right(20u32)
        .margin_bottom(10u32)
        .margin_left(10u32)
        .x_label_area_size(40u32)
        .y_label_area_size(60u32)
        .build_cartesian_2d(0f64..x_max, 0f64..y_max)?;

    chart
        .configure_mesh()
        .bold_line_style(TRANSPARENT)
        .light_line_style(TRANSPARENT)
        .x_desc("Putative reads (M)")
        .y_desc("Distinct fragments (M)")
        .axis_desc_style(("sans-serif", 12u32))
        .label_style(("sans-serif", 10u32))
        .draw()?;

    chart.draw_series(LineSeries::new(data.clone(), GREEN.stroke_width(2)))?;
    chart.draw_series(
        data.iter()
            .map(|&pt| Circle::new(pt, 3u32, GREEN.filled())),
    )?;

    plot_area.present()?;
    Ok(())
}
