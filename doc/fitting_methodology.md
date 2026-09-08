# Spectrum-fitting methodology of `scripts/` — formalism and scientific context

_Disassembly of `scripts/{config,gasmodel,instrument,fitting,dataio}.py` (v0.1.0), written 2026-09-02.
Every equation below was checked against the source and, where numerically checkable, against a
re-derivation run on the actual package — see §9._

---

## 1. What the code is doing, in one sentence

The package inverts continuum-subtracted SPHEREx spectrophotometry of a cometary coma for the three
production rates $Q_{\mathrm{H_2O}}$, $Q_{\mathrm{CO_2}}$, $Q_{\mathrm{CO}}$ by writing the coma
emission as a **strictly linear** function of those rates,

$$F_\nu(\lambda_i) \;=\; \sum_X Q_X\, A_{iX},$$

and solving the resulting weighted linear least-squares problem analytically. Everything upstream of
that solve exists to build the design matrix $A$; everything downstream exists to decide what part
of the solution may honestly be quoted.

The scientific lineage is explicit in `config.py`: this is the AKARI/IRC method of
**Ootsubo et al. (2012), ApJ 752, 15** — simultaneous 2.7 / 4.3 / 4.7 µm band photometry of a coma,
converted to production rates through fluorescence efficiencies from Crovisier's molecular database
(**Bockelée-Morvan & Crovisier 1989**; band model of **Crovisier & Encrenaz 1983**) — transplanted
onto SPHEREx's linear-variable-filter spectrophotometry, and turned from a per-band conversion into
a single joint multi-species fit.

---

## 2. Layer architecture

The design is a strict pipeline; each layer is a separate, individually testable function.

| Layer | Module | Content | Key symbol |
|---|---|---|---|
| 0 | `gasmodel` | geometry, $r_h$-scaled lifetimes, expansion velocity | $\tau_X(r_h)$, $v_g(r_h)$ |
| 1 | (implicit) | Haser parent density | $n_X(r)$ |
| 2 | `gasmodel` | aperture integration, Yamamoto filling factor | $f_1(x)$, $N_{\rm ap}$ |
| 3 | `gasmodel` | fluorescence excitation, optional pump opacity | $g_b(r_h)$ |
| 4 | `gasmodel` | monochromatic emission spectrum | $\mathcal{S}_X(\lambda)$ |
| 5 | `instrument` | SPHEREx channel bandpass | $W_{ij}$ |
| 6 | `fitting` | design matrix, WLS solve, errors, limits | $A_{iX}$, $\hat{Q}$, $C$ |

`dataio` sits outside the physics: it assembles the fit input from `data/emission/` and persists
`results/`.

---

## 3. Layers 0–2: geometry to column density

### 3.1 Lifetimes and expansion velocity

Photodissociation lifetimes are scaled with the inverse-square dilution of the solar UV field,

$$\tau_X(r_h) \;=\; \tau_X(1\,\mathrm{au})\; r_h^{2},$$

with `TAU_1AU = {H2O: 8.3e4, CO2: 5.0e5, CO: 1.3e6}` s (`config.TAU_1AU`, Ootsubo et al. 2012
Table 2). The expansion velocity is **fixed, not fitted**, to the empirical law

$$v_g(r_h) \;=\; 0.8\, r_h^{-1/2}\ \mathrm{km\,s^{-1}}$$

(`ModelParams.v_g`; overridable with `v_g_kms`). [CITATION NEEDED: source of the
$0.8\,r_h^{-1/2}$ expansion-velocity law — the constant is hard-coded and unattributed in
`config.py`.]

The docstring gives the reason $v_g$ is frozen: in the small-aperture regime the band flux goes as
$Q/v_g$, so the two are *exactly* degenerate and a joint fit is unbounded. §7.2 quantifies what that
choice costs.

### 3.2 Haser parent coma

A single-scale-length (parent-only) Haser profile,

$$n_X(r) \;=\; \frac{Q_X}{4\pi v_g r^{2}}\,\exp\!\left(-\frac{r}{v_g \tau_X}\right),$$

whose line-of-sight column at impact parameter $\rho$ in the thin-coma limit is
$N_X(\rho) = Q_X/(4 v_g \rho)$. This form appears *explicitly* in the code only inside
`opacity_factor` mode 2; in the default optically thin path it enters only through its aperture
integral, §3.3. There are no daughter species and no dissociation-product bookkeeping — the model
describes the parent molecules alone, which is the appropriate level of description for the
vibrational fundamentals being observed.

### 3.3 Aperture integration — the Yamamoto filling factor

For a centred circular aperture of projected radius $\rho_{\rm ap}$ the number of parent molecules
in the beam is

$$\boxed{\;N_{\rm ap} \;=\; Q_X\,\tau_X\, f_1(x), \qquad x \equiv \frac{\rho_{\rm ap}}{v_g \tau_X}\;}$$

with the Yamamoto (1981) filling factor

$$f_1(x) \;=\; x\,g(x), \qquad
g(x) \;=\; \frac{1}{x} - K_1(x) + \int_x^{\infty} K_0(y)\,\mathrm{d}y ,$$

$K_0, K_1$ being modified Bessel functions of the second kind. `filling_factor()` implements this in
three numerical regimes, which is the only interesting thing about the function:

- **$x < 10^{-4}$** — the term $1/x - K_1(x)$ is a catastrophic cancellation of two $\mathcal{O}(1/x)$
  quantities, so the series is used instead:
  $$g(x) = \frac{\pi}{2} + x\left[\tfrac12\ln\frac{x}{2} + \tfrac{\gamma}{2} - \tfrac34\right] + \mathcal{O}(x^3),$$
  $\gamma$ the Euler–Mascheroni constant.
- **$10^{-4} \le x \le 30$** — closed form via modified Struve functions,
  $$\int_x^\infty K_0(y)\,\mathrm{d}y = \frac{\pi}{2}\Big[1 - x\big(K_0(x)\mathbf{L}_{-1}(x) + K_1(x)\mathbf{L}_0(x)\big)\Big].$$
- **$x > 30$** — asymptotic series
  $\sqrt{\pi/2x}\,e^{-x}\left(1 - \tfrac{5}{8x} + \tfrac{129}{128x^2}\right)$, because the closed
  form loses precision to cancellation there.

`method="quad"` re-does the integral by direct quadrature as an independent cross-check; the two
agree to $1.2\times10^{-13}$ relative (§9).

**The two limits carry all the physics.**

$$f_1(x) \to \frac{\pi x}{2}\ (x\to0) \qquad\Longrightarrow\qquad N_{\rm ap} \to \frac{\pi\,Q_X\,\rho_{\rm ap}}{2\,v_g},$$

the classic small-aperture result — **independent of the lifetime**, because the aperture samples
only the innermost coma where no molecule has yet had time to dissociate. And $f_1 \to 1$
($x\to\infty$), where the aperture contains the entire coma and $N_{\rm ap} = Q\tau$.

Which regime a given species is in matters a great deal, and the three species are *not* in the same
one. At $\rho_{\rm ap} = 40\,000$ km:

| $r_h$ | $x_{\rm H_2O}$ | $f_1$ | $x_{\rm CO_2}$ | $f_1$ | $x_{\rm CO}$ | $f_1$ |
|---|---|---|---|---|---|---|
| 1.0 au | 0.602 | 0.555 | 0.100 | 0.138 | 0.039 | 0.057 |
| 1.5 au | 0.328 | 0.368 | 0.054 | 0.079 | 0.021 | 0.032 |
| 2.5 au | 0.152 | 0.199 | 0.025 | 0.038 | 0.010 | 0.015 |
| 4.0 au | 0.075 | 0.106 | 0.013 | 0.019 | 0.005 | 0.008 |

CO and CO₂ are always deep in the small-aperture limit; H₂O, with its 60× shorter lifetime, sits in
the *transition* region at typical SPHEREx distances. This asymmetry propagates into the mixing
ratios — see §7.2.

### 3.4 The geometric amplitude

`amplitude_per_unit_Q` packages Layers 0–2 plus the $r_h^{-2}$ of the excitation and the
$1/4\pi\Delta^2$ of the inverse square law into a single scalar per species per channel:

$$\mathcal{A}_X(r_h, \Delta) \;=\; \frac{\tau_X(r_h)\, f_1(x_X)}{4\pi \Delta^2}\; r_h^{-2},
\qquad \Delta\ \text{in metres}.$$

Note that **each channel keeps its own $r_h$ and $\Delta$** (`build_design_matrix`). SPHEREx scans a
moving target non-simultaneously, so within one "phase" the geometry varies by a few per cent;
collapsing to a mean would bias the amplitude by the same amount.

---

## 4. Layer 3: excitation

### 4.1 Fluorescence efficiencies

Band $g$-factors scale with the solar flux only:

$$g_b(r_h) \;=\; \frac{g_b(1\,\mathrm{au})}{r_h^{2}} \qquad [\mathrm{photons\;s^{-1}\;molecule^{-1}}]$$

Eight bands are tabulated (`config.BANDS`, Ootsubo et al. 2012 Table 2):

| Species | Band | $\lambda_b$ (µm) | $g_b(1\,\mathrm{au})$ | hot |
|---|---|---|---|---|
| H₂O | $\nu_3$ | 2.66 | 2.82e-4 | |
| H₂O | $\nu_1$ | 2.73 | 2.47e-5 | |
| H₂O | $\nu_2+\nu_3-\nu_2$ | 2.66 | 2.81e-5 | ✓ |
| H₂O | $\nu_1+\nu_3-\nu_1$ | 2.73 | 2.19e-5 | ✓ |
| H₂O | $\nu_3-\nu_2$ | 4.63 | 7.66e-6 | ✓ |
| H₂O | $\nu_1-\nu_2$ | 4.85 | 7.35e-6 | ✓ |
| CO₂ | $\nu_3$ | 4.26 | 2.9e-3 | |
| CO | $v(1\!-\!0)$ | 4.67 | 2.6e-4 | |

This is a **pure resonance-fluorescence, optically thin** excitation model: no collisional
excitation, no dependence of $g$ on rotational or vibrational temperature, no radiation-trapping
correction. `T_rot` (default 70 K) is accepted throughout but currently only reaches the band
*profile*, not the band *strength*.

### 4.2 Optional pump-opacity corrections

`opacity_factor(species, Q, r_h, params)` returns a multiplicative correction $\eta_X$. Three modes:

- **Mode 0 (default)** — $\eta = 1$ exactly. This is the *only* mode that leaves the model linear in
  $Q$, and therefore the only one usable inside the fit.
- **Mode 1 — core excision.** An opaque core of radius $\rho_\tau = \rho_{\tau,\rm ref}(Q/10^{27})$
  is removed from the aperture:
  $$\eta = \frac{f_1(x_{\rm ap}) - f_1(x_\tau)}{f_1(x_{\rm ap})}, \qquad x_\tau = \frac{\rho_\tau}{v_g\tau}.$$
  $\rho_{\tau,\rm ref} = 10^3$ km is a **placeholder** calibrated for CO₂ $\nu_3$ at $r_h=1.3$ au
  (Debout et al. 2016, Icarus 265, 110).
- **Mode 2 — first-order pump correction** in the PSG manner (Villanueva et al.). With the pump
  optical depth $\tau_p(\rho) = \kappa_p N(\rho)$ built from the thin-coma column and a
  Sobolev-like escape probability $\beta(\tau_p) = (1-e^{-\tau_p})/\tau_p$,
  $$\eta = \frac{1}{\rho_{\rm ap}-\rho_0}\int_{\rho_0}^{\rho_{\rm ap}} \beta\big(\tau_p(\rho)\big)\,\mathrm{d}\rho.$$
  $\kappa_p \equiv S_p/v_p$ is a **placeholder** anchored so $\tau_p \sim 1$ at $\rho \sim 1000$ km
  for $Q(\mathrm{CO_2}) = 10^{27}\,$s$^{-1}$, with the other bands scaled by $g$.

Both non-default modes are documented as diagnostics for *mapping* the systematic after a fit, never
for use inside one. See §8 for a caveat on the mode-2 average.

---

## 5. Layer 4: the monochromatic spectrum, and why the fit is linear

Each band's total flux at Earth is the standard optically thin emission integral,

$$F_b \;=\; \frac{N_{\rm ap}\; g_b(r_h)\; (hc/\lambda_b)}{4\pi \Delta^2} \qquad [\mathrm{W\,m^{-2}}],$$

distributed over a **unit-area** profile $\Phi_b(\lambda)$ ($\int\Phi_b\,\mathrm{d}\lambda = 1$):

$$F_\lambda(\lambda) \;=\; \sum_b F_b\, \Phi_b(\lambda), \qquad
F_\nu[\mathrm{mJy}] \;=\; F_\lambda[\mathrm{W\,m^{-2}\,\mu m^{-1}}]\; \lambda_{\mu m}^2 \times \frac{10^{23}}{c}.$$

`band_profile()` is currently a Gaussian of the band's tabulated intrinsic FWHM — a **placeholder**
for real unit-normalised fluorescence templates from PSG / the GSFC Fluorescence Database or
Crovisier's LESIA line lists. Because the profiles are unit-normalised, replacing them changes band
*shapes* but leaves every band-integrated flux — and hence every retrieved $Q$ — unchanged at
SPHEREx resolution.

**The structural fact that makes everything downstream linear.** For one species all bands share the
same $N_{\rm ap}$, the same $r_h^{-2}$ and the same $1/4\pi\Delta^2$. The *relative* weights of a
species' bands are therefore fixed regardless of geometry, and the spectrum factorises exactly:

$$\boxed{\;F_{\nu,X}(\lambda) \;=\; Q_X \cdot \underbrace{\mathcal{A}_X(r_h,\Delta)}_{\text{geometry}} \cdot \underbrace{\mathcal{S}_X(\lambda)}_{\text{fixed shape}}\;}$$

$$\mathcal{S}_X(\lambda) \;=\; \left[\sum_{b\in X} g_b(1\,\mathrm{au})\,\frac{hc}{\lambda_b}\,\Phi_b(\lambda)\right] \lambda_{\mu m}^2 \frac{10^{23}}{c}\quad[\mathrm{mJy\ per\ unit\ amplitude}].$$

Each species contributes one geometry-independent template times one scalar. That is precisely the
structure of a linear model with three free parameters.

---

## 6. Layer 5: the instrument

SPHEREx disperses with linear variable filters, so every measurement is a flux density through one
filter channel, characterised in the photometry table by its centre `wl` and FWHM `wlwidth`. The
module builds a normalised Gaussian throughput

$$T_i(\lambda) = \frac{1}{\sigma_i\sqrt{2\pi}} \exp\left[-\frac{(\lambda-\lambda_i)^2}{2\sigma_i^2}\right],
\qquad \sigma_i = \frac{\mathrm{FWHM}_i}{2\sqrt{2\ln 2}},$$

and band-averages the model exactly as the detector does:

$$\bar{F}_i \;=\; \frac{\int F_\nu(\lambda)\,T_i(\lambda)\,\mathrm{d}\lambda}{\int T_i(\lambda)\,\mathrm{d}\lambda}
\;\;\longrightarrow\;\; \bar{F}_i = \sum_j W_{ij}\,F_\nu(\lambda_j), \qquad \sum_j W_{ij} = 1.$$

Two implementation choices worth noting. The rows of $W$ are normalised with the **trapezoid weights
of the internal grid itself**, so the matrix product is a proper weighted mean under the same
integration rule used everywhere else; and the profile is truncated at $5\sigma$ before
renormalisation, which omits $\sim6\times10^{-7}$ of the weight without biasing the mean. Averaging
is done in $F_\nu$ rather than $F_\lambda$, which avoids ever having to define an effective
wavelength for the conversion.

The internal grid is log-spaced at constant resolving power $R = 4000$ over 2.0–5.2 µm — roughly 30×
finer than the instrument, so the band-average is not grid-limited.

Measured from the photometry itself, the SPHEREx channels in use have

| Region | median $\lambda/\Delta\lambda$ |
|---|---|
| 2.3 µm | ≈ 40 |
| 2.4–3.0 µm (H₂O window) | ≈ 35 |
| 4.0–5.0 µm (CO₂ / CO windows) | ≈ 124 |

At $R\approx35$ the 2.7 µm complex is a single unresolved blend of four H₂O bands; at $R\approx124$
the 4.26 µm CO₂ $\nu_3$, the 4.67 µm CO fundamental and the 4.63 / 4.85 µm H₂O hot bands are
marginally separated at best. This is the whole reason the fit must be joint (§7.3), and the reason
the concept document flags the as-built LSF as the highest-priority placeholder.

---

## 7. Layer 6: the inverse problem

### 7.1 Design matrix and solve

For channel $i$ and species $X$,

$$A_{iX} \;=\; \mathcal{A}_X(r_{h,i}, \Delta_i)\;\sum_j W_{ij}\,\mathcal{S}_X(\lambda_j)
\qquad [\mathrm{mJy\ per\ molecule\ s^{-1}}],$$

i.e. the flux one molecule per second of species $X$ would deposit in that channel, at that
channel's own geometry and through its own filter. The observation model is then

$$d_i \;=\; \sum_X A_{iX} Q_X + \varepsilon_i, \qquad \varepsilon_i \sim \mathcal{N}(0, \sigma_i^2),$$

with $d_i$ = `emis_mjy` and $\sigma_i$ = `emis_err_mjy`, and the estimator is the ordinary weighted
least-squares minimiser of

$$\chi^2(\mathbf{Q}) \;=\; \sum_i \frac{1}{\sigma_i^2}\left(d_i - \sum_X A_{iX}Q_X\right)^2 .$$

Numerically the system is solved on the whitened, column-scaled matrix
$\tilde{A}_{iX} = (A_{iX}/\sigma_i)/s_X$ with $s_X = \max_i |A_{iX}/\sigma_i|$, because raw entries
are $\sim10^{-28}$ mJy per molecule s$^{-1}$ while $Q \sim 10^{28}$ — an unscaled normal matrix would
be numerically hopeless. The covariance is un-scaled afterwards:

$$C \;=\; \frac{(\tilde{A}^{\mathsf T}\tilde{A})^{+}}{s\,s^{\mathsf T}}\, \cdot\, \kappa^2,
\qquad \kappa = \begin{cases}\sqrt{\chi^2_\nu} & \chi^2_\nu > 1 \ \text{and rescaling enabled}\\ 1 & \text{otherwise}\end{cases}$$

with $\nu = n_{\rm chan} - n_{\rm free}$. The $\chi^2_\nu$ rescaling is the standard remedy for
formal photometric errors that understate the true scatter — here by up to $\sim60\times$, because
one "phase" stitches together scans taken days apart. `Q_err_formal` preserves the unrescaled value
so the two can be compared. Note the rescaling is **one-sided**: $\chi^2_\nu < 1$ never shrinks the
errors, which is the conservative choice.

**The solve is unconstrained.** A species consistent with zero can and does come back negative; that
is reported rather than clipped, because forcing $Q \ge 0$ would bias every non-detection upward.

### 7.2 Sensitivity to the fixed parameters

Because $f_1$ is locally a power law, $f_1 \propto x^{s}$ with $s \equiv \mathrm{d}\ln f_1/\mathrm{d}\ln x$,
the retrieved rate for a fixed measured flux scales as

$$Q_X \;\propto\; v_g^{\,s_X}\; \tau_X^{\,s_X - 1}\; \rho_{\rm ap}^{-s_X}.$$

$s \to 1$ in the small-aperture limit and $s \to 0$ when the coma is fully enclosed. Evaluated at
$\rho_{\rm ap} = 40\,000$ km, $r_h = 1$ au: $s_{\rm H_2O} = 0.61$, $s_{\rm CO_2} = 0.89$,
$s_{\rm CO} = 0.95$; by $r_h = 2.5$ au all three have risen towards the small-aperture limit
(0.85 / 0.96 / 0.98).

Three consequences that belong in any error budget:

1. **$v_g$ is the dominant systematic in absolute $Q$**, as the docstring says — but not quite at
   full strength for water: at $r_h = 1$ au a 20 % error in $v_g$ moves $Q(\mathrm{CO_2})$ by
   17.7 % and $Q(\mathrm{H_2O})$ by only 11.7 %.
2. **It does not cancel in the mixing ratio**, contrary to the usual assumption. Because the species
   sit in different regimes,
   $$\frac{Q_{\rm CO_2}}{Q_{\rm H_2O}} \propto v_g^{\,s_{\rm CO_2}-s_{\rm H_2O}} \approx v_g^{0.29}\ (r_h = 1\,\mathrm{au}),$$
   so a 20 % $v_g$ error leaves a residual $+5.4$ % in the ratio, falling to $+2.0$ % at
   $r_h = 2.5$ au. Small, but real, and heliocentric-distance dependent: it is a *shrinkage* of the
   systematic rather than a cancellation, and it is worth stating explicitly rather than assuming
   zero.
3. **Lifetime errors touch water and only water.** $Q \propto \tau^{s-1}$ vanishes as $s\to1$, so
   $Q(\mathrm{CO})$ and $Q(\mathrm{CO_2})$ are essentially lifetime-independent at these apertures,
   while $Q(\mathrm{H_2O}) \propto \tau_{\rm H_2O}^{-0.39}$ at 1 au. A 30 % lifetime error is a
   $-9.9$ % water systematic — and, since $\tau \propto r_h^2$, it introduces a residual $r_h$
   dependence into the H₂O/CO₂ ratio that the other two species do not share.

### 7.3 Why the fit must be joint

The 4.55–4.90 µm window contains CO $v(1\!-\!0)$ at 4.67 µm blended with the H₂O hot bands
$\nu_3-\nu_2$ (4.63 µm) and $\nu_1-\nu_2$ (4.85 µm). Taking the tabulated $g$-factors at face value
and a typical CO/H₂O mixing ratio of $\sim1$ %, the water hot bands are several times brighter than
the CO fundamental in that window. Fitting the region in isolation would overestimate
$Q(\mathrm{CO})$ by close to an order of magnitude. Fitting it *simultaneously* with 2.7 µm — where
$Q(\mathrm{H_2O})$ is pinned by the strong $\nu_3$/$\nu_1$ complex, which then predicts the hot-band
contribution at 4.6–4.9 µm — is what makes $Q(\mathrm{CO})$ meaningful at all.

The code encodes this as the `h2o_anchored` flag: `True` iff the 2.7 µm band entered the fit. The
docstring is blunt that this flag, not the error bar, is what tells you whether to believe
$Q(\mathrm{CO})$.

### 7.4 Coverage, effective channel count, and detection classification

Three separate guards sit between the raw solve and a quotable number.

**(a) Coverage.** Coverage is decided on narrow *diagnostic* ranges (`KEY_RANGES`), tighter than the
continuum-subtraction windows: H₂O 2.60–2.80 µm (≥3 channels), CO₂ 4.20–4.30 µm (≥2), CO
4.60–4.70 µm (≥2). The distinction being drawn is categorical and is the methodologically strongest
idea in the module:

> *not covered* $\neq$ *non-detection*. A non-detection means the band was observed and no flux was
> found; not-covered means the band was never adequately sampled, so **no statement about $Q$ — not
> even an upper limit — is supportable.**

Counts are taken on the channels that survive filtering, since an excluded channel constrains
nothing. A species can have a perfectly healthy model column from the wings of a neighbouring band
while its own diagnostic feature was never sampled; the key-range test catches exactly that case.

**(b) Effective channel count.** The participation ratio of the per-channel Fisher information,

$$w_{iX} = \left(\frac{A_{iX}}{\sigma_i}\right)^{2}, \qquad
n_{\rm eff,X} = \frac{\left(\sum_i w_{iX}\right)^{2}}{\sum_i w_{iX}^{2}},$$

which equals $n$ for $n$ equally-informative channels and $\to 1$ when one channel dominates. This
is the guard against formally significant values carried by a single measurement — the handoff notes
a 24σ $Q(\mathrm{H_2O})$ of $10^{28}$ s$^{-1}$ at $r_h = 3.4$ au resting on one channel at the edge
of the 2.7 µm window. The rule of thumb encoded in `caveats()` is $n_{\rm eff} \ge 2$.

**(c) Classification.** With $k$ = `upper_limit_sigma`:

| Condition | `status` | reported |
|---|---|---|
| $n_{\rm key} <$ `min_points` | `not_covered` | nothing |
| $\hat{Q} < 0$ | `negative_fit` | limit $= k\sigma_Q$ only |
| $0 \le \hat{Q} < k\sigma_Q$ | `upper_limit` | $\hat{Q}$ and limit $= \hat{Q} + k\sigma_Q$ |
| $\hat{Q} \ge k\sigma_Q$ | `detected` | $\hat{Q} \pm \sigma_Q$ |

The `negative_fit` limit is the $\max(\hat{Q},0) + k\sigma$ convention with $\hat{Q}<0$.

**Channel selection** (`channel_masks`, single source of truth so the figures mark exactly what the
solve dropped): require finite $d_i$, finite $\sigma_i > 0$; drop $d_i < -k_-\sigma_i$ with
$k_- = 1$; optionally apply a symmetric $|d_i| > k_c\sigma_i$ clip (disabled by default). The
negative cut is deliberately one-sided — a significantly negative flux is an over-subtracted
continuum, not a measurement of negative emission, whereas a positive outlier may be real.

### 7.5 Mixing ratios

$$r = \frac{Q_a}{Q_b}, \qquad
\sigma_r^2 = r^2\left(\frac{C_{aa}}{Q_a^2} + \frac{C_{bb}}{Q_b^2} - \frac{2C_{ab}}{Q_a Q_b}\right),$$

i.e. full first-order propagation *including* the covariance term — which matters, because
$Q(\mathrm{H_2O})$ and $Q(\mathrm{CO})$ are strongly anti-correlated through the 4.6–4.9 µm blend.
Ratios are computed only when both species are `detected`; a negative best fit has no production
rate and cannot enter a ratio.

---

## 8. Upstream: the continuum subtraction

The fitter's input, `emis_mjy`, is produced by `notebooks/continuum_subtraction.ipynb`, not by this
package. The model there is a local weighted polynomial in $\lambda-\lambda_c$ per band,

$$c(\lambda) = \sum_{m=0}^{p} a_m (\lambda-\lambda_c)^m, \qquad
\sigma_c(\lambda) = \sqrt{\mathbf{v}(\lambda)^{\mathsf T} C_a \mathbf{v}(\lambda)},$$

$\mathbf{v}$ the Vandermonde row, fitted on continuum windows with the emission windows punched out,
with:

- $p$ requested per band, **capped at 1 whenever the continuum sample does not bracket the emission
  window** — curvature you have to extrapolate is curvature you cannot justify;
- a 2-point fallback (nearest blue + nearest red) when too few continuum points survive;
- leave-one-out (or 10-fold) cross-validated order selection as a *check*, reported as
  `cv_best_order`;
- validation on reduced $\chi^2$, `err_scale` $=\mathrm{rms}/\mathrm{median}(\sigma)$, a two-halves
  residual-shape $z$-score, CV tolerance, bracketing, and the significance of any negative excursion
  of the continuum across the emission window;
- $d = F - c$, $\sigma_d = \sqrt{\sigma_F^2 + \sigma_c^2}$.

Only `role == "emission"` channels from bands with an accepted verdict (`PASS`/`WARN`) enter the fit.
A `FAIL` band's residual is the residual of a continuum the validation rejected, so it is not a
measurement of anything.

Physically, the continuum being removed is the sum of scattered sunlight from dust and thermal
emission from the same grains, which cross over in the 3–4 µm region; a low-order local polynomial
is a deliberately agnostic stand-in for a dust model (cf. the Harker / Kelley references in
`doc/dust-grain/`). The price is that any broad *spectral* feature — the 3 µm water-ice band in
particular (cf. `doc/water-ice/`, Protopapa et al. 2014) — is partly absorbed into the continuum and
partly into the residual.

---

## 9. Verification performed

Run against the installed package on the real catalogue, 2026-09-02:

| Check | Result |
|---|---|
| $f_1(x) \to \pi x/2$ as $x\to0$ | ratio 0.999995 at $x=10^{-6}$ |
| $f_1(x) \to 1$ as $x\to\infty$ | 1.0000 at $x=10^{5}$ |
| Struve closed form vs. direct quadrature | max relative difference $1.2\times10^{-13}$ |
| Band-integrated $F_\nu$ from `spectrum_mjy`, converted back to W m⁻², vs. the closed-form $N_{\rm ap}\sum_b g_b(hc/\lambda_b)/4\pi\Delta^2$ | agree to $1\times10^{-8}$ |
| `WM2UM_TO_MJY` $=10^{23}/c$ | dimensionally correct for $F_\nu[\mathrm{mJy}] = F_\lambda\lambda_{\mu m}^2\cdot10^{23}/c$ |
| Linearity in $Q$ | $F(3Q)/F(Q) = 3.0000000$ |
| Bandpass row normalisation | $\sum_j W_{ij} = 1$ exactly; $W\mathbf{1} = \mathbf{1}$ |
| WLS solution vs. independent normal-equation solve | identical to 9 significant figures ($\hat{Q}$ and $\sigma_{\rm formal}$) |
| $n_{\rm eff}$ vs. independent participation ratio | identical |

---

## 10. Assumptions, limitations and issues found

### 10.1 Physical assumptions (by construction)

1. **Pure resonance fluorescence, optically thin.** No collisional excitation anywhere; $g$ scales as
   $r_h^{-2}$ with no temperature dependence. The hot-band $g$-factors in particular depend on the
   vibrational excitation of the inner coma, which is *not* a pure $r_h^{-2}$ quantity — and the
   CO/H₂O separation at 4.7 µm rests entirely on those hot-band strengths. This is a physics
   limitation of the tabulated-$g$ approach, not a coding gap.
2. **Single-scale-length Haser, parent only**, spherically symmetric, steady state, isotropic
   outflow, constant $v_g$ with radius. No jets, no extended sources, no nucleus.
3. **Constant $v_g$ across species.** All three molecules are given the same expansion velocity.
4. **Aperture centred on the nucleus**, circular, and $\rho_{\rm ap}$ exactly matching the
   photometry — enforced only by convention (`ModelParams.rho_ap_km` "must match the `r_ap_km` of
   the input photometry"), not by an assertion.
5. **Gaussian LSF and Gaussian band envelopes**, both flagged placeholders.

### 10.2 Statistical issues found in the code

**(a) The continuum uncertainty is correlated but treated as diagonal.** `emis_err_mjy` folds
$\sigma_c(\lambda)$ into the per-channel error in quadrature, and the fit then uses a diagonal weight
matrix. But $\sigma_c$ comes from a *single* polynomial fit per band, so its errors are strongly
correlated across the channels of that band — the notebook's own `_trapz_band_flux` acknowledges
this and sums the continuum term **linearly** rather than in quadrature. The fitter does not. The
consequence is that $\sigma_Q$ underestimates the continuum-systematic contribution, in the
direction where it matters most: a coherent tilt or offset of the continuum under a band is exactly
degenerate with a change in $Q$.

*This is fixable with what is already on disk.* The emission CSVs carry `cont_c0..c3` and the
coefficient covariance, so the data covariance can be built as

$$C_{d,ij} = \sigma_{F,i}^2\,\delta_{ij} + \delta_{b(i)b(j)}\;\mathbf{v}_i^{\mathsf T} C_a^{(b)} \mathbf{v}_j,$$

block-diagonal by band, and the solve becomes the generalised least squares
$\hat{Q} = (A^{\mathsf T} C_d^{-1} A)^{-1} A^{\mathsf T} C_d^{-1} d$ with
$C_Q = (A^{\mathsf T}C_d^{-1}A)^{-1}$. Recommended as the single highest-value change to §7.

**(b) `upper_limit_sigma` defaults to 1.0, so `status == "detected"` means "≥1σ".** In a live fit
run during this review, a species with $\hat{Q}/\sigma_Q = 1.38$ was labelled `detected`. Nothing is
wrong internally — the number and its error are both reported — but `detected` is a load-bearing
word: it gates entry into `ratio()` and will be read by anyone using the catalogue as meaning a
detection. Either raise the default to 3, or rename the tier.

**(c) The one-sided negative-channel cut biases upward** (already in `handoff.md`, +0.288σ per
channel for pure noise). It interacts with (b): the same channels that push a marginal $\hat{Q}$
above $1\sigma$ are the ones the cut preferentially retains. Worth quantifying jointly with a
noise-only injection test before quoting any 1–3σ rate.

**(d) `clip_sigma` is a signal cut, not an outlier cut.** `channel_masks` drops channels with
$|d_i| > k_c\sigma_i$ — computed on the data, not on the fit residual. Enabling it would remove the
*brightest genuine emission channels* first. It is `None` by default, and should probably be
re-specified as a residual cut (a second pass on $d - A\hat{Q}$) before anyone reaches for it.

**(e) Opacity mode 2 averages the escape probability linearly in $\rho$.**
`np.trapezoid(corr, rr)/(rr[-1]-rr[0])` weights each impact parameter equally, whereas an aperture
average weights by $2\pi\rho\,\mathrm{d}\rho$. Since $\beta$ is smallest (most opaque) at small
$\rho$, the linear average over-weights the opaque core and overestimates the correction. Additive
fix: weight by `rr` in both integrals. The grid is also linear from 1 km to $\rho_{\rm ap}$, so the
inner coma where $\tau_p$ actually matters is sampled by a handful of points out of 2000 — a
log-spaced grid would be better. Mode 2 is diagnostic-only, so this is not currently affecting any
published number.

**(f) `dof` counts only covered species.** $\nu = n_{\rm chan} - n_{\rm free}$ with
$n_{\rm free} = \sum_X \mathbb{1}[\text{covered}]$, which is correct for the solve — but it means
$\chi^2_\nu$, and therefore the error rescaling $\kappa$, is computed on a model that may be missing
a species whose flux is physically present in the data. The `caveats()` machinery flags the specific
case that matters (H₂O dropped while CO is fitted, so the 4.63 µm hot band goes unmodelled into
$Q(\mathrm{CO})$), but the inflated $\chi^2_\nu$ then *widens* every error bar rather than being
attributed to the missing component. Interpret $\kappa \gg 1$ as a model-incompleteness flag, not
just as noise.

### 10.3 Placeholders blocking publication-grade numbers

From `config.PLACEHOLDERS` and the handoff, in priority order:

1. **The SPHEREx instrument model** — as-built LVF centres, $R(\lambda)$, and the real LSF. The
   CO/H₂O separability at 4.7 µm depends directly on $R$ there, and the Gaussian is only guaranteed
   to have the right *width*, not the right wings.
2. **Band profiles** $\Phi_b$ — real unit-area fluorescence templates at $T_{\rm rot}$. Neutral for
   $Q$ by the unit-normalisation argument, but required before any claim about band *shape*.
3. `RHO_TAU_REF_KM` and `KAPPA_PUMP` — needed only to map the opacity systematic, not to fit.
4. The $v_g$ law's provenance (§3.1).

---

## 11. Summary

The methodology is a clean, deliberately linear implementation of the Ootsubo et al. (2012)
band-fluorescence method, with two genuine methodological advances over a per-band conversion: the
**joint three-species solve**, which is the only way to break the CO / H₂O-hot-band blend at
4.6–4.9 µm, and the **explicit separation of "not covered" from "not detected"**, which prevents the
catalogue from reporting limits that the data cannot support. The numerics are careful — three
asymptotic regimes for $f_1$, an independent quadrature cross-check, column scaling before the normal
equations, per-channel geometry — and every equation checked here reproduces to machine precision.

The two things to fix before the catalogue is quoted are statistical rather than physical: the
correlated continuum error currently propagated as if diagonal (§10.2a), and the 1σ definition of
`detected` (§10.2b). The two things that will limit the *accuracy* rather than the precision are
$v_g$ (§7.2) and the as-built instrument model (§10.3).
