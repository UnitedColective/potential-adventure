# PCF Methodology & Technical Notes

## Mathematical Foundation

### The Second Difference

For any time series X, we define:
- **First difference**: ΔXₜ = Xₜ₊₁ − Xₜ
- **Second difference**: Δ²Xₜ = ΔXₜ₊₁ − ΔXₜ = Xₜ₊₂ − 2Xₜ₊₁ + Xₜ

The second difference is a natural discrete approximation of the second derivative (acceleration) and is commonly used to remove linear trends while preserving higher-order structure.

### Lag-1 Autocorrelation

For a series Y with mean μ and variance σ²:

ρ₁(Y) = E[(Yₜ − μ)(Yₜ₊₁ − μ)] / σ²

In practice, we estimate using:

ρ̂₁(Y) = Σ(Yₜ − Ȳ)(Yₜ₊₁ − Ȳ) / Σ(Yₜ − Ȳ)²

### White-Noise Baseline

For pure white noise Nₜ ~ N(0, 1):

**Key result**: E[ρ₁(Δ²N)] = −2/3 exactly

This can be derived algebraically:
- Δ²Nₜ = Nₜ₊₂ − 2Nₜ₊₁ + Nₜ
- Var(Δ²N) = E[(Δ²N)²] = 6 (since E[Nᵢ²] = 1 and E[NᵢNⱼ] = 0 for i ≠ j)
- Cov(Δ²Nₜ, Δ²Nₜ₊₁) involves products of 4 Gaussian terms
- After calculation: ρ₁(Δ²N) = −2/3

### The PCF Metric

δ = ρ₁(Δ²X) + 2/3

Properties:
- **δ = 0** under white noise (by construction)
- **δ > 0** indicates positive autocorrelation in second differences (memory/structure)
- **δ < 0** indicates negative autocorrelation in second differences (over-damping)
- **|δ|** scales with intensity of temporal structure

## Statistical Testing

### Significance Against White Noise

Given an observed δ, we test the null hypothesis H₀: X ~ white noise.

**Method**: Generate n white noise samples of length N, compute δ for each, and compute:

z = (δ_obs − μ_wn) / σ_wn

where μ_wn and σ_wn are the empirical mean and SD under white noise.

**P-value**: Two-tailed test using standard normal.

### Surrogate Testing

Two standard null models:

#### 1. Shuffle Surrogate
- Randomly permute the time series
- Destroys temporal order; retains distribution
- Tests for temporal coherence

#### 2. Phase-Randomized Surrogate
- Compute FFT of the original series
- Randomize phase while preserving amplitude spectrum
- Preserves spectral content; destroys phase coupling

For each surrogate method:
- Generate n surrogates (current: n=100; target: n≥1000)
- Compute δ for each surrogate
- Compute z-score: z = (δ_obs − mean(δ_surr)) / std(δ_surr)

## Confidence Intervals

### Block Bootstrap

We use **stationary block bootstrap** to account for temporal dependence:

1. Choose block size b = ⌊√N⌋
2. Randomly select blocks of length b from the series
3. Concatenate blocks until length N is reached
4. Compute δ on resampled series
5. Repeat m times (typically m = 1000)
6. Estimate CI via percentiles

**Why block bootstrap?** Standard bootstrap assumes IID data, which is violated in time series.

## Domain-Specific Notes

### Prime Gaps

**Data**: Gap between consecutive primes: gₙ = pₙ₊₁ − pₙ

**Challenges**:
- Distribution is non-stationary (gaps grow with magnitude of primes)
- Must check stability across sample sizes
- Residue structure is discrete (mod 24)

**Checks performed**:
- Stability: verified δ consistent for N ∈ {10⁵, 10⁶, 3×10⁸}
- Detrending: δ changes by <0.001 after removing low-frequency trend

### Cardiac RR Intervals

**Data**: Beat-to-beat interval series from ECG

**Challenges**:
- Highly dependent on heart rate, respiration, autonomic state
- Medication effects may bias results
- Within-group variance is high

**Analysis design**:
- Extract first 5,000 valid RR intervals (minimizes drift)
- Test for stratification by heart rate
- Check robustness with windowed analysis

**Expected pattern**:
- Healthy: δ > 0 (adaptive variability)
- Heart failure: δ ≈ 0 (loss of complexity) — *not strongly confirmed*

### Neural iEEG

**Data**: Intracranial recording from amygdala during fear-conditioning task

**State classification**: Assign each window to E (excitable), H (hybrid), or B (bistable) based on eigenspectrum.

**Challenges**:
- Small sample (n=3 subjects, 17 trials each)
- Different metric used (composite stability score, not δ directly)
- Requires electrode localization and preprocessing

**Result**: Superlinear relationship S ~ f^1.746 suggests hybrid state is an attractor that strengthens nonlinearly with occupancy.

### Circadian Gene Expression

**Data**: Microarray probe intensity over 8 circadian timepoints

**Challenges**:
- Only 8 timepoints (insufficient for robust δ estimation)
- High dimensionality (44,544 probes) requires dimension reduction
- Must account for batch effects

**Approach**: Fit Fourier model to each probe; compare phase alignment distribution between MCF7 (cancer) and MCF10A (normal).

**Result**: Cancer cells show *higher* phase coherence (counterintuitive), suggesting rigidification rather than loss of rhythm.

## Robustness Checks

All analyses include:

1. **Detrending**: Remove polynomial trend; recompute δ
2. **Winsorizing**: Trim outliers at 1% and 99%; recompute δ
3. **Windowed analysis**: Compute δ over sliding windows; check stability
4. **Surrogate testing**: Compare to shuffled and phase-randomized nulls

## Caveats and Limitations

1. **Sample sizes**: Cardiac (n=54, 29), Neural (n=3), Circadian (n=44,544 probes but 8 timepoints)
2. **Surrogate count**: Current 100, should be ≥1,000 per standard practice
3. **Multiple comparisons**: No correction applied (exploratory analysis)
4. **Causality**: PCF is correlational; cannot infer mechanistic causation
5. **Generalization**: Each domain uses different preprocessing; direct comparison limited

## Future Methodological Work

- [ ] Increase surrogate count to 1,000–10,000
- [ ] Develop domain-specific preprocessing standards
- [ ] Test for scale-invariance (multifractal analysis)
- [ ] Extend to multivariate time series
- [ ] Derive theoretical connection to information theory / criticality
- [ ] Develop clinical decision thresholds (e.g., for CHF risk)

---

*For implementation details, see `code/pcf_core.py`*
