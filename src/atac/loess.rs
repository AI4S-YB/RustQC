//! Minimal loess port matching stats::loess.smooth defaults
//! (span=2/3, degree=2, family="gaussian", no robust reweighting).

/// Fit local degree-2 polynomial weighted by tricube on the q nearest neighbors,
/// where q = ceil(span * n). Evaluates at every x in `xs`.
pub fn loess_smooth(xs: &[f64], ys: &[f64], span: f64, degree: usize) -> Vec<f64> {
    assert_eq!(xs.len(), ys.len());
    let n = xs.len();
    if n == 0 {
        return vec![];
    }
    let q = ((span * n as f64).ceil() as usize).clamp(degree + 1, n);
    let mut out = Vec::with_capacity(n);
    for &x0 in xs {
        // Pick q nearest neighbors by |x - x0|.
        let mut dists: Vec<(usize, f64)> = (0..n).map(|i| (i, (xs[i] - x0).abs())).collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let nbrs: Vec<usize> = dists.iter().take(q).map(|(i, _)| *i).collect();
        let max_d = dists[q - 1].1.max(f64::MIN_POSITIVE);
        // Tricube weights.
        let w: Vec<f64> = nbrs
            .iter()
            .map(|&i| {
                let u = (xs[i] - x0).abs() / max_d;
                let one_minus = (1.0 - u.powi(3)).max(0.0);
                one_minus.powi(3)
            })
            .collect();
        // Solve weighted least squares y ~ poly(x − x0, degree) by normal equations.
        // Build X (q × (degree+1)) and W (diagonal, weights).
        let p = degree + 1;
        let mut xtwx = vec![0.0f64; p * p];
        let mut xtwy = vec![0.0f64; p];
        for (k, &i) in nbrs.iter().enumerate() {
            let dx = xs[i] - x0;
            let mut row = vec![1.0f64; p];
            for j in 1..p {
                row[j] = row[j - 1] * dx;
            }
            let wk = w[k];
            for a in 0..p {
                for b in 0..p {
                    xtwx[a * p + b] += row[a] * row[b] * wk;
                }
                xtwy[a] += row[a] * ys[i] * wk;
            }
        }
        // Solve (p × p) symmetric positive (semi-)definite system via Gauss-Jordan.
        let beta = solve_linear(&mut xtwx, &mut xtwy, p);
        // Fitted value at x0 corresponds to the constant term β₀.
        out.push(beta[0]);
    }
    out
}

fn solve_linear(a: &mut [f64], b: &mut [f64], p: usize) -> Vec<f64> {
    // In-place Gaussian elimination on (a | b).
    for k in 0..p {
        // Pivot.
        let mut piv = k;
        for r in k + 1..p {
            if a[r * p + k].abs() > a[piv * p + k].abs() {
                piv = r;
            }
        }
        if piv != k {
            for c in 0..p {
                a.swap(k * p + c, piv * p + c);
            }
            b.swap(k, piv);
        }
        let akk = a[k * p + k];
        if akk.abs() < 1e-15 {
            return vec![0.0; p];
        }
        for r in 0..p {
            if r == k {
                continue;
            }
            let factor = a[r * p + k] / akk;
            for c in k..p {
                a[r * p + c] -= factor * a[k * p + c];
            }
            b[r] -= factor * b[k];
        }
    }
    (0..p).map(|i| b[i] / a[i * p + i]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_quadratic_exactly_at_full_span() {
        let xs: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 + 3.0 * x + 0.5 * x * x).collect();
        // span=1.0 + degree=2 → exact recovery of a quadratic.
        let fit = loess_smooth(&xs, &ys, 1.0, 2);
        for (a, b) in fit.iter().zip(ys.iter()) {
            assert!(
                (a - b).abs() < 1e-6,
                "loess(span=1) on quadratic: {} vs {}",
                a,
                b
            );
        }
    }

    #[test]
    fn fits_constant_signal() {
        let xs: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let ys = vec![5.0; 20];
        let fit = loess_smooth(&xs, &ys, 2.0 / 3.0, 2);
        for v in fit {
            assert!((v - 5.0).abs() < 1e-9, "{}", v);
        }
    }
}
