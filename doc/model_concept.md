# Synthetic Coma Gas Emission Model for SPHEREx (0.7–5.0 µm)

**Concept design document**
Scope: $\mathrm{H_2O}$, $\mathrm{CO_2}$, $\mathrm{CO}$ fluorescence emission; low spectral resolution ($R \approx 40$–$100$); forward model from gas parameters to observed spectrum.
Source literature: `doc/gas-fluorescence/` — Crovisier & Encrenaz (1983), Ootsubo et al. (2012), Kelley et al. (2016), Debout et al. (2016), Villanueva et al. (2011, 2018), PSG Handbook.

* * *

## 1\. Fundamental theories

### 1\.1 Resonance fluorescence in a collisionless medium

The coma is optically thin and collisionless over almost all of its volume. Molecules absorb solar photons and re\-emit before any collision redistributes the energy, so level populations are set by **fluorescence equilibrium**, not by a local kinetic temperature. This is the regime that makes the whole model tractable: the emission per molecule depends only on the solar radiation field, not on the local density.

The emission rate per molecule in vibrational band *b* is the **fluorescence efficiency**, or band $g$-factor:

$$g_b(r_h) = g_b(1\text{ au}) \cdot r_h^{-2} \quad \left[\text{photons s}^{-1}\text{ molecule}^{-1}\right]$$

The $r_h^{-2}$ scaling is simply solar flux dilution. Physically, $g_b$ is the sum over all pump transitions of the solar photon absorption rate multiplied by the branching ratio into band $b$:

$$g_b = \sum_{\text{pumps}} \frac{F_\odot(\nu_{lu})}{h\nu_{lu}} \cdot \sigma_{lu} \cdot B(lu \to b)$$

Crovisier & Encrenaz (1983) computed these for cometary parent molecules using a vibrational band model; Villanueva et al. (2011, 2012) extended the treatment to full line\-by\-line cascade, which matters for $\mathrm{H_2O}$ and $\mathrm{NH_3}$ where cascade contributions are large. **Two second\-order dependences exist and can be neglected at $R \approx 40$–$100$:** the rotational temperature (which redistributes intensity among lines within a band but preserves the band total) and the heliocentric velocity (Swings effect, which shifts pump lines against Fraunhofer structure). Neither changes a band\-integrated flux by more than a few percent.

### 1\.2 Haser coma structure

Steady, spherically symmetric outflow at constant velocity with exponential photodissociation:

$$n_X(\rho) = \frac{Q_X}{4\pi v_g \rho^2} \exp\!\left[-\frac{\rho}{v_g \tau_X}\right]$$

This is Ootsubo et al. (2012) Eq. 1 exactly. The photodissociation lifetime scales with solar flux:

$$\tau_X(r_h) = \tau_X(1\text{ au}) \cdot r_h^2$$

and the expansion velocity follows the standard empirical law

$$v_g = 0.8 \cdot r_h^{-0.5} \quad [\text{km s}^{-1}]$$

though $v_g$ is a **free parameter** in this model (see §4) precisely because it is the largest single systematic in absolute $Q$.

### 1\.3 Column density and aperture integration

The observable is not the local density but the number of molecules within the beam. Integrating the Haser density along the line of sight and over a circular aperture of radius $\rho_{ap}$ gives the number of molecules in the field of view:

$$N_{ap} = Q_X \cdot \tau_X \cdot f_1(x), \qquad x = \frac{\rho_{ap}}{v_g \tau_X}$$

where $f_1$ is the **cometary filling factor** for a parent species in a centred cylindrical beam (Yamamoto 1981, as implemented in PSG):

$$f_1(x) = x \cdot g(x)$$

$$g(x) = \frac{1}{x} - K_1(x) + \int_x^{\infty} K_0(y)\, dy$$

with $K_0$, $K_1$ the modified Bessel functions of the second kind. The two limits are the sanity checks worth hard\-coding as unit tests:

- **Small aperture** ($x \to 0$, $\rho_{ap} \ll v_g \tau$): $g(x) \to \pi/2$, so $N_{ap} \to \dfrac{\pi Q \rho_{ap}}{2 v_g}$. Photodissociation is irrelevant; the column density is the pure $1/\rho$ profile.
- **Large aperture** ($x \to \infty$): $f_1 \to 1$, $N_{ap} \to Q\tau$. The beam contains the entire coma.

For SPHEREx at typical comet distances, apertures of $10^3$–$10^5$ km sit in the small-$x$ regime for $\mathrm{CO_2}$ ($\tau = 5.0\times10^5$ s) and CO ($\tau = 1.3\times10^6$ s), and closer to intermediate for $\mathrm{H_2O}$ ($\tau = 8.3\times10^4$ s). The full Bessel form should be used throughout rather than the $1/\rho$ approximation — $\mathrm{H_2O}$ is the species where the difference bites.

### 1\.4 From column to flux

The band flux received at the observer:

$$F_b = N_{ap} \cdot g_b(r_h) \cdot \frac{hc}{\lambda_b} \cdot \frac{1}{4\pi\Delta^2} \quad [\text{W m}^{-2}]$$

In the small\-aperture limit this reduces to Kelley et al. (2016) Eq. 4:

$$F_b = \frac{Q\, \rho_{ap}\, h c\, g_b}{8\, \lambda_b\, v_g\, \Delta^2\, r_h^2}$$

The two forms are algebraically identical once $N_{ap} = \pi Q \rho_{ap}/(2v_g)$ is substituted, which is a useful cross\-check on your normalization. Verify the factor of $\pi$ explicitly during implementation — this is the single most common place to lose a 3.14 in this class of model.

### 1\.5 Optical depth: the Ootsubo assumption and its correction

Ootsubo et al. (2012) adopt the **optically thin approximation**, justified not by claiming the coma is everywhere thin but by the argument that *"the inner coma region, which is optically thick, is small and does not significantly contribute to the observed flux of the coma."* This holds when $\rho_{ap} \gg \rho_\tau$, where $\rho_\tau$ is the radius inside which the solar pump is attenuated. For AKARI's $1' \times 1'$ aperture ($10^4$–$10^5$ km at the comet) against an opaque core of $10^1$–$10^3$ km, the assumption is safe.

Debout et al. (2016) quantify $\rho_\tau$ for the $\mathrm{CO_2}\ \nu_3$ band at $r_h = 1.3$ au through full 3\-D radiative transfer:

| $Q(\mathrm{CO_2})\ [\text{s}^{-1}]$ | $g_{\text{eff}}/g_{\text{thin}}$ at nucleus | Converges to $g_{\text{thin}}$ at |
| --- | --- | --- |
| $10^{25}$ | $\approx 1.00$ (thin everywhere) | — |
| $10^{26}$ | $\approx 0.90$ | $\rho \approx 100$ km |
| $10^{27}$ | $\approx 0.30$ | $\rho \approx 1000$ km |

So $\rho_\tau$ scales roughly linearly with $Q$, and the opaque core stays below \~$10^3$ km even for a very active comet. This is the numerical justification for adopting Ootsubo's assumption as the model default.

PSG implements a **first\-order analytic correction** worth building in as an optional mode. For each $g$-factor it tracks a weighted representative pump line intensity $S_p$; the pump optical depth across the column is

$$\tau_p = N_{col} \cdot \frac{S_p}{v_p}$$

and the effective efficiency becomes

$$g_{\text{thick}} = g_{\text{thin}} \cdot \frac{1 - \exp(-\tau_p)}{\tau_p}$$

valid at low solar phase angle. This has the right asymptotics ($\to g_{\text{thin}}$ for $\tau_p \to 0$, $\to g_{\text{thin}}/\tau_p$ for $\tau_p \gg 1$) and costs almost nothing to evaluate.

**Design decision:** implement opacity as a switchable module with three modes — (0) fully thin, the Ootsubo default; (1) core excision, zeroing emission inside a radius $\rho_\tau(Q)$; (2) the PSG multiplicative correction applied radially. Mode 0 is the production default for the catalog; modes 1–2 quantify the systematic.

### 1\.6 Instrumental convolution

At $R = 40$–$100$ no ro\-vibrational line is resolved. Each band collapses toward a single feature whose width is set by the instrument, not by the molecule. But band *widths* are not all negligible: the $\mathrm{H_2O}$ 2.7 µm complex ($\nu_1 + \nu_3$ plus hot bands) spans roughly 0.1 µm, comparable to the resolution element at $R \approx 35$–$40$, while $\mathrm{CO_2}\ \nu_3$ at 4.26 µm is intrinsically narrow and will appear as an unresolved delta function. The model must therefore carry a band *profile*, not just a band total.

* * *

## 2\. Model concepts

### 2\.1 Architecture

A six\-layer forward pipeline, each layer independently testable:

**Layer 0 — Geometry.** From ($r_h$, $\Delta$, aperture in km), or from an ephemeris query. Produces the projected aperture radius $\rho_{ap}$ and the solar flux scaling.

**Layer 1 — Coma density.** Haser profile per species from ($Q_X$, $v_g$, $\tau_X(r_h)$).

**Layer 2 — Aperture integration.** $N_{ap,X} = Q_X \tau_X f_1(x_X)$ via the Bessel filling factor. This is where the aperture dependence enters, and it is *species\-dependent* because $\tau$ differs by more than an order of magnitude between $\mathrm{H_2O}$ and CO.

**Layer 3 — Excitation.** $g_b(r_h) = g_b(1\text{ au})/r_h^2$, optionally corrected for pump opacity. Returns an effective per\-band photon emission rate.

**Layer 4 — Monochromatic emission.** Build the high\-resolution radiance, $F_\lambda^{\text{gas}}(\lambda)$, by distributing each band's total flux over a normalized band profile. Two implementations, discussed in §2.2.

**Layer 5 — Instrument.** Convolve with the line spread function at the local $R(\lambda)$, then resample onto the SPHEREx wavelength grid.

**Layer 6 — Output.** $F_\lambda(\lambda)$ over 0.7–5.0 µm, with per\-band integrated fluxes retained as diagnostics.

Layers 0–3 are analytic and cheap. Layer 4 carries all the physics complexity. Layers 5–6 are instrument bookkeeping.

### 2\.2 Band profiles: the key implementation choice

**Option A — Pure band model (Ootsubo\-style).** Represent each band as a delta function at $\lambda_b$ carrying flux $F_b$, then convolve. Fastest, and formally adequate for integrated flux, but it produces the wrong band shape for $\mathrm{H_2O}$ 2.7 µm and cannot represent the hot\-band structure near 4.6–4.9 µm.

**Option B — Line\-by\-line (PSG/GSFC).** Generate the full fluorescence line list at a chosen $T_{\text{rot}}$, then convolve. Physically correct, but expensive to run inside a fitting loop over thousands of catalog objects.

**Recommended: Option C — precomputed normalized templates.** Run PSG (or the GSFC fluorescence database) once to generate high\-resolution spectra for each species on a grid of ($T_{\text{rot}}$, $r_h$), normalize each to unit total band emission, and store as templates. The runtime model then becomes a fast linear combination:

$$F_\lambda^{\text{gas}}(\lambda) = \sum_X \sum_b \left[ N_{ap,X} \cdot g_b(r_h) \cdot \frac{hc}{\lambda_b} \cdot \frac{1}{4\pi\Delta^2} \right] \cdot \Phi_b(\lambda; T_{\text{rot}})$$

where $\Phi_b$ is the stored unit\-normalized profile. This preserves band shape and hot\-band structure at the cost of one precomputation, and keeps the runtime cost of a synthetic spectrum at the level of a few array operations. It also makes the model trivially invertible, since $F$ is linear in each $Q$.

### 2\.3 The dominant low\-resolution problem: $\mathrm{H_2O}$ hot\-band contamination of CO

This is not a refinement. It is the first\-order structural issue for a SPHEREx gas model, and it must be built into the forward model rather than applied as a post\-hoc correction.

The CO $v(1\text{–}0)$ fundamental at 4.67 µm sits between two $\mathrm{H_2O}$ hot bands: $\nu_3-\nu_2$ at 4.63 µm and $\nu_1-\nu_2$ at 4.85 µm. At $R \approx 100$ these are marginally separable; at $R \approx 40$ they are not. Using Ootsubo's Table 2 $g$-factors, consider a comet with a typical $\mathrm{CO}/\mathrm{H_2O} = 1\%$:

- CO contribution: $Q_{CO} \cdot 2.6\times10^{-4} = 2.6\times10^{-6} \cdot Q_{H_2O}$
- $\mathrm{H_2O}$ hot bands: $Q_{H_2O} \cdot (7.66+7.35)\times10^{-6} = 1.5\times10^{-5} \cdot Q_{H_2O}$

**The water hot bands are roughly six times brighter than the CO fundamental itself.** Ootsubo et al. handle this by subtracting the hot\-band contribution using the ratio of hot\-band to 2.7 µm $g$-factors, which is only possible because they measure the 2.7 µm band in the same spectrum. Your model must do the same, but as a joint forward model: the 4.6–4.9 µm region is a *blend* whose $\mathrm{H_2O}$ component is tied to $Q(\mathrm{H_2O})$ retrieved at 2.7 µm. Treating the 4.67 µm feature as pure CO will overestimate $Q(\mathrm{CO})$ by nearly an order of magnitude at low mixing ratios.

The direct consequence for the catalog: **$Q(\mathrm{CO})$ is not independently measurable from the 4.7 µm region alone.** It is a differential measurement against a well\-constrained $Q(\mathrm{H_2O})$, and its uncertainty inherits everything that goes wrong at 2.7 µm — continuum placement, ice\-band wing, and $\mathrm{H_2O}$ opacity.

### 2\.4 Structure of the 0.7–5.0 µm range

With only $\mathrm{H_2O}$, $\mathrm{CO_2}$ and CO, the gas component of the model is **identically zero below \~2.5 µm**. The overtone and combination bands of water near 1.38 and 1.9 µm have $g$-factors orders of magnitude below the fundamentals and are undetectable in fluorescence. This is a feature, not a limitation: the 0.7–2.5 µm region is pure dust continuum, and it provides the lever arm needed to extrapolate the continuum underneath the 2.7 µm water band — which is exactly the hardest continuum placement in the whole range.

The emitting range therefore partitions as:

| Region | Content |
| --- | --- |
| 0\.7 – 2.5 µm | Dust scattering continuum only; continuum anchor |
| 2\.5 – 2.9 µm | $\mathrm{H_2O}$ $\nu_3+\nu_1+$ hot bands; sits on the 3 µm ice\-band wing |
| 4\.2 – 4.3 µm | $\mathrm{CO_2}\ \nu_3$; narrow, strong, first to go optically thick |
| 4\.6 – 4.9 µm | CO $v(1\text{–}0)$ blended with $\mathrm{H_2O}\ \nu_3-\nu_2$ and $\nu_1-\nu_2$ |

Note the asymmetry in instrumental resolution across these features: the $\mathrm{H_2O}$ band falls in SPHEREx's lowest-$R$ bands while $\mathrm{CO_2}$ and CO fall in the higher-$R$ long\-wavelength bands. Confirm the exact $R(\lambda)$ breakpoints against the instrument specification before finalizing the LSF model — this asymmetry works in your favour for $\mathrm{CO}/\mathrm{CO_2}$ and against you for $\mathrm{H_2O}$.

* * *

## 3\. Assumptions

Stated explicitly, each with its validity limit:

**A1. Spherically symmetric, isotropic outgassing.** Real comae show sunward enhancement (a factor of a few in 3I/ATLAS). Acceptable for aperture\-integrated fluxes; fails for any spatially resolved analysis.

**A2. Constant expansion velocity, independent of radius and species.** Gas actually accelerates through the inner coma and heavy/light species differ. The error is largely absorbed into the fitted $v_g$.

**A3. Single\-generation Haser decay with exponential lifetime.** Ootsubo et al. estimate the resulting systematic at \~10% for large apertures and include it in their error budget. Adopt the same 10% floor.

**A4. Fluorescence equilibrium everywhere in the beam.** Breaks down inside the collisional radius, \~200 km for $Q(\mathrm{H_2O}) = 10^{28}\text{ s}^{-1}$, scaling roughly linearly with $Q$. Negligible for apertures $\gg 10^3$ km.

**A5. Optically thin emission** (see §1.5). The default. Quantified and switchable.

**A6. No extended (distributed) sources.** Sublimating icy grains produce $\mathrm{H_2O}$ far from the nucleus, and both 3I/ATLAS and C/2017 K2 show a non\-asymptoting $Q(\mathrm{H_2O})$ curve. In a fixed\-aperture, spatially unresolved model this appears as an *apparent* $Q(\mathrm{H_2O})$ that grows with aperture size. Document it as a known bias rather than attempting to model it in v1.

**A7. Fixed rotational temperature.** $T_{\text{rot}}$ affects intra\-band distribution only, which is unresolved at $R \leq 100$. Adopt a fixed value (50–80 K is conventional) and record it. Verify by generating templates at $T_{\text{rot}} = 30, 50, 100$ K and confirming the convolved band shapes agree within noise.

**A8. Fluorescence equilibrium $g$-factors scale exactly as $r_h^{-2}$.** Ignores the Swings effect. Sub\-percent for band totals.

**A9. Nucleus and dust contribute no gas emission**, and the gas model adds linearly to the dust continuum model. Standard, and correct as long as the coma is thin.

* * *

## 4\. Parameters and output

### 4\.1 Free parameters

| Symbol | Quantity | Units | Typical range |
| --- | --- | --- | --- |
| $Q_{H_2O}$, $Q_{CO_2}$, $Q_{CO}$ | Gas production rates | $\text{molecules s}^{-1}$ | $10^{23}$ – $10^{30}$ |
| $r_h$ | Heliocentric distance | au | 0\.5 – 10 |
| $v_g$ | Gas expansion velocity | $\text{km s}^{-1}$ | 0\.2 – 1.0 |
| $\rho_{ap}$ | Aperture radius at the comet | km | $10^2$ – $10^6$ |
| $R$ | Spectral resolving power | — | 40 – 130 |

Note that $\Delta$ (observer distance) is needed to convert $\rho_{ap}$ from angular to physical units and to compute the $1/\Delta^2$ dilution. Either treat $\Delta$ as an additional geometric input or supply $\rho_{ap}$ directly in km and $\Delta$ separately — do not conflate them.

### 4\.2 Fixed (dependent) parameters

Band $g$-factors at 1 au and photodissociation lifetimes, following Ootsubo et al. (2012) Table 2, taken from Bockelée\-Morvan & Crovisier (1989) and Crovisier's molecular database:

| Molecule | Lifetime (s) | Band | $\lambda$ (µm) | $g$-factor ($\text{s}^{-1}$) |
| --- | --- | --- | --- | --- |
| $\mathrm{H_2O}$ | $8.3\times10^4$ | $\nu_3$ | 2\.66 | $2.82\times10^{-4}$ |
|  |  | $\nu_1$ | 2\.73 | $2.47\times10^{-5}$ |
|  |  | $\nu_2+\nu_3-\nu_2$ | 2\.66 | $2.81\times10^{-5}$ |
|  |  | $\nu_1+\nu_3-\nu_1$ | 2\.73 | $2.19\times10^{-5}$ |
|  |  | $\nu_3-\nu_2$ | 4\.63 | $7.66\times10^{-6}$ |
|  |  | $\nu_1-\nu_2$ | 4\.85 | $7.35\times10^{-6}$ |
| $\mathrm{CO_2}$ | $5.0\times10^5$ | $\nu_3$ | 4\.26 | $2.9\times10^{-3}$ |
| CO | $1.3\times10^6$ | $v(1\text{–}0)$ | 4\.67 | $2.6\times10^{-4}$ |

Lifetimes are quoted at 1 au and scale as $r_h^2$. Also fixed: $T_{\text{rot}}$, the band profiles $\Phi_b$, and the solar\-pump reference spectrum implicit in the $g$-factors.

### 4\.3 Intermediate diagnostics to expose

Worth returning alongside the spectrum, because they are what you will need when a fit misbehaves: $x_X$ and $f_1(x_X)$ per species; $N_{ap,X}$; $\tau_p$ per band; $g_{\text{eff}}/g_{\text{thin}}$ per band; and the per\-band integrated flux $F_b$ before convolution.

### 4\.4 Output

Primary: $F_\lambda(\lambda)$ [$\text{W m}^{-2}\,\mu\text{m}^{-1}$] or $F_\nu$ [Jy] on the SPHEREx wavelength grid across 0.7–5.0 µm, gas component only, to be added to the dust continuum and thermal model.

Because $F$ is strictly linear in each $Q$ at fixed ($r_h$, $v_g$, $\rho_{ap}$), the natural output structure is a **design matrix**: three unit-$Q$ template spectra that the fitting stage combines. This makes retrieval a linear least\-squares problem for the $Q$'s at fixed geometry, with only the continuum requiring nonlinear treatment.

* * *

## 5\. Required input data

### 5\.1 Have (in `gas-fluorescence/`)

- **Ootsubo et al. 2012, ApJ 752, 15** — Table 2 $g$-factors and lifetimes; the Haser formalism; the optically\-thin justification; Table 3 provides 18 comets of validation data at exactly this wavelength range and comparable resolution.
- **Kelley et al. 2016, PASP 128, 018009** — the band\-flux equation and the two\-component coma spectrum framework.
- **Debout et al. 2016, Icarus 265, 110** — the quantitative opacity calibration in §1.5.
- **Villanueva et al. 2011, Icarus 216 / 2018, JQSRT 217, 86 / PSG Handbook** — fluorescence model construction, the Yamamoto filling factor, the $g_{\text{thick}}$ correction.
- **Crovisier & Encrenaz 1983, A&A 126, 170** — the original band $g$-factor computation. *Note: this PDF is a scanned image with no extractable text layer; OCR it or work from the tabulated values reproduced in Ootsubo Table 2.*

### 5\.2 Still needed

**Fluorescence line lists / band profiles.** The GSFC Fluorescence Database via PSG, or Crovisier's LESIA database (`lesia.obspm.fr/perso/jacques-crovisier/`). Required for Option C templates. PSG can be driven through its API to generate the templates in batch.

**Aperture integration code.** PSG's `isocoma` package implements the Yamamoto filling factors, including offset beams and square pixels: `github.com/nasapsg/isocoma`. Worth adopting rather than reimplementing the Bessel integrals.

**Photodissociation rates.** Huebner & Mukherjee (2015) for active vs. quiet Sun, if you want lifetimes beyond Ootsubo's single set. The active/quiet difference is a factor of \~2 in $\tau$ and propagates directly into $N_{ap}$ in the large-$x$ regime.

**SPHEREx instrument model.** The LVF band centres, $R(\lambda)$ per band, and throughput. This is the single largest gap — everything in Layer 5 depends on it, and the $\mathrm{H_2O}$-versus-$\mathrm{CO_2}/\mathrm{CO}$ resolution asymmetry noted in §2.4 cannot be assessed without it.

**Solar reference spectrum.** Only if recomputing $g$-factors from scratch rather than adopting tabulated values; Kurucz plus ACE, as used by PSG.

### 5\.3 Validation targets

Ootsubo et al. Table 3 is the natural benchmark: 18 comets with $Q(\mathrm{H_2O})$, $Q(\mathrm{CO_2})$, $Q(\mathrm{CO})$ derived from 2.5–5 µm spectra at low resolution using exactly this formalism. Reproducing their production rates from their measured band fluxes is a closed\-loop test of Layers 1–4. Lisse et al. (2025, 2026) provide the only published SPHEREx comet spectra (3I/ATLAS) and are the end\-to\-end test including Layer 5.

* * *

## 6\. Recommended build order

1. Layers 0–3 with the band\-model approximation and Ootsubo's $g$-factors. Validate against Ootsubo Table 3 — you should recover their $Q$ values from their fluxes to within their quoted errors.
2. Add the Bessel filling factor and confirm the two analytic limits in §1.3 as unit tests.
3. Generate Option C band\-profile templates from PSG; confirm that convolving them at $R = 40$ reproduces the band\-model fluxes.
4. Implement the joint $\mathrm{H_2O}/\mathrm{CO}$ treatment at 4.6–4.9 µm (§2.3) and quantify the $Q(\mathrm{CO})$ bias from ignoring it.
5. Add the opacity module and map $g_{\text{eff}}/g_{\text{thin}}$ across the ($Q$, $\rho_{ap}$, $r_h$) parameter space to define where the Ootsubo assumption is safe for the catalog.
6. Couple to the dust continuum and thermal model.
