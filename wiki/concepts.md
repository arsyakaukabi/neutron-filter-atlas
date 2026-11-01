# Concepts

What this page is for: the definitions a reader needs to interpret the numbers of the atlas, with their formulas.
Who should read it: every user of the atlas, before comparing crystals.

The values behind every definition (bands, spectra, thresholds) are in `configs/`, mostly
`configs/figures_of_merit.yaml` and `configs/uncertainty.yaml`; [Configuration](configuration.md) lists them. Symbols:
λ wavelength, E energy, n atom number density, xᵢ atomic fraction of element i, σ cross section per atom,
μ = n·σ attenuation coefficient (atoms/Å³ × barn = 1/cm).

## Thermal band and spectrum

A filter is judged by what it does to the thermal neutrons the instrument wants. The atlas weights wavelengths in the
**thermal band** 0.5–3.0 Å by a Maxwellian flux from a moderator at 300 K:

φ(λ) ∝ λ⁻⁵ exp(−λ_T² / λ²),  λ_T = h / √(2 mₙ k_B · 300 K) = 1.77889 Å.

The **thermal transmission** of a slab of thickness L is the flux-weighted average over the band (200-point
trapezoid rule):

T_th(L) = ∫ φ(λ) exp(−μ(λ) L) dλ / ∫ φ(λ) dλ.

A second, **hot band** of 0.3–1.0 Å gives `L08_hot`, with the same 300 K Maxwellian weights (not a hot-source
spectrum).

μ(λ) is the **off-Bragg** attenuation of a large single crystal oriented so that no Bragg reflection is excited:
inelastic scattering (phonons), incoherent elastic scattering and absorption, computed by NCrystal. Coherent elastic
(Bragg) scattering is left out. Real crystals show narrow Bragg dips on top of this baseline at particular
wavelengths and orientations; see [Comparing crystals](comparing-crystals.md#what-the-table-cannot-tell).

## Filter thickness L08

**L₀.₈** (column `L08`) is the thickness at which the off-Bragg crystal transmits 80 % of the thermal beam:
T_th(L₀.₈) = 0.80. It is found by a root search for every crystal and temperature.

Why a fixed transmission: a filter has to be thin enough to keep the thermal beam and thick enough to remove the
rest. Fixing the thermal cost (every crystal loses the same 20 % of the band) turns the comparison into one question:
how much unwanted radiation does each crystal remove at the thickness where it costs exactly that much? A crystal
with low thermal attenuation can be made thick, which is why L₀.₈ itself is not the figure of merit but the
thickness at which the figures of merit are evaluated.

## Rejections

The column names call every neutron rejection `A_fast`, the epithermal one included. The figures of merit are
**rejections**: A = −ln T of an unwanted component at L₀.₈. A rejection of 1 means a factor
e, 2.3 a factor 10, 4.6 a factor 100. Each comes with `_lo` and `_hi`, the same rejection at L₀.₈·exp(∓σ_model)
(see [Model uncertainty](#model-uncertainty)).

| Column | Component | Definition |
|---|---|---|
| `A_fast_1E` | epithermal neutrons | 1/E spectrum, 1 eV–100 keV (the slowing-down spectrum of a moderator) |
| `A_fast_watt` | fission neutrons | Watt spectrum S(E) ∝ exp(−E/a) sinh(√(bE)), a = 0.988 MeV, b = 2.249 /MeV, over 100 keV–10 MeV |
| `A_gamma_1MeV` | 1 MeV photons | μ_γ(1 MeV) · L₀.₈ |
| `A_gamma_6MeV` | 6 MeV photons | μ_γ(6 MeV) · L₀.₈ |

For neutrons: A_fast(L) = −ln[ ∫ S(E) exp(−n σ_tot(E) L) dE / ∫ S(E) dE ], with σ_tot = Σ xᵢ σ_tot,i the free-atom
total cross section from ENDF/B-VIII.0 at 294 K (natural isotopic abundances; a grid of 5,000 points per decade).
For photons: μ_γ = n Σ xᵢ σ_γ,i with coherent and incoherent scattering, photoelectric absorption and pair
production.

**Uncollided, not all transmitted.** These rejections count only particles that cross the slab without any
interaction (narrow-beam attenuation). A real filter also passes some particles that scatter and still leave the
back face, so the rejection of all transmitted particles is lower, by an amount that depends on the beam
collimation and geometry. For the shortlist the pipeline computes the all-particle value with OpenMC slab transport
(`stage4_transport.json`). Read the atlas rejections as the values for a well-collimated beam; they rank filters
optimistically for an uncollimated beam.

Two single-energy summaries are also given: `mu_epi_1eV`, the free-atom attenuation at 1 eV, and `nieman_ratio` =
`mu_epi_1eV / mu_1.8`, the epithermal attenuation per unit of thermal attenuation at 1.8 Å.

## Three levels of phonon description

The thermal attenuation depends on the lattice dynamics, because phonons scatter thermal neutrons inelastically. The
atlas computes it at three levels of detail:

1. **Composition proxy** (`proxy_L08`, `proxy_mu_1.8`, `proxy_T_1cm`): no phonons at all.
   μ_proxy(λ) = n Σ xᵢ [σ_abs,i · λ/1.798 Å + σ_inc,i], free-element absorption (1/v law) and incoherent scattering.
   It leaves out inelastic scattering, so it is a lower bound on μ and `proxy_L08` an upper bound on L₀.₈ (up to a
   small tolerance measured on the crystals with both values). `proxy_L08` uses the density of the Materials Project
   summary structure; `L08_proxy_vdos_n` is the same quantity at the density of the full model, for like-for-like
   comparison with `L08`.
2. **Elastic Debye** (`debye_L08`, `theta_D_elastic`): one Debye temperature for every atom of the crystal,
   θ_D from the elastic moduli of the Materials Project (Anderson's formula with the Voigt-Reuss-Hill average). Only
   crystals with an elastic tensor have it.
3. **Full phonon model** (`L08` and every rejection): per-element vibrational densities of states from the phonon
   density of states, either the Materials Project DFT phonons (`vdos_source == "mp_dft"`) or phonons computed for
   the atlas with the learned interatomic potential SevenNet (`"sevennet"`). DOS values below 2 meV are dropped and
   NCrystal extrapolates g(E) ∝ E² below that energy.

Every crystal has level 1. Level 2 exists where `debye_L08` is not empty. Level 3 exists where `vdos_source` is not
`none`; only level 3 carries rejections, flags of the lattice and a model uncertainty.

## Coverage: which crystals have a full phonon model

`funnel_step` records how far each crystal went:

| `funnel_step` | Meaning |
|---|---|
| `mp_dft` | full model from Materials Project DFT phonons |
| `sevennet` | full model from SevenNet phonons |
| `excluded_by_proxy_bound` | even the composition-only upper bound transmits less than half of the band through 1 cm, so the crystal cannot be a filter |
| `unstable_dft` | the Materials Project phonons are dynamically unstable (imaginary modes), so they cannot be used |
| `mp_dos_conversion_failed` | the Materials Project DOS could not be converted to per-element DOS |
| `not_computed_budget` | no Materials Project phonons, and outside the fixed budget of 1,500 learned-potential calculations (chosen in the order of `priority`) |
| `supercell_too_large` | the phonon supercell would exceed 480 atoms |
| `mlip_failed` | the learned-potential relaxation did not converge within 500 steps or the phonon calculation failed |

Crystals without a full model keep their proxy and, where available, Debye values;
`stage3_not_computed_table.csv` lists them with the position their Debye value or their upper bound would give
them ([Results tables](results-tables.md#crystals-without-a-full-model)).

## Tiers: commercially grown crystals and other candidates

`tier == "A"` marks the **commercially grown crystals**: entries matched to a vendor table of crystals offered for
sale (`data/curated/growable_crystals.csv`; how each was matched is in `tierA_match`). `tier == "B"` is every other
candidate. A tier-B crystal may be excellent on paper and still unavailable as a large single crystal.

The candidate set itself: experimentally known Materials Project entries within 0.05 eV/atom of the convex hull,
without hydrogen, radioactive elements or strong thermal absorbers (B, Cd, Gd, Sm, Eu, Li, Hg, In, Ir, Rh, Au, Hf, Dy,
Er, Co, Ag); see `configs/candidates.yaml`.

## Model uncertainty

The model is compared with measured total cross-sections (see [Validation](#validation-against-measured-cross-sections)),
and the scatter of that comparison becomes a model uncertainty on ln L₀.₈, `sigma_model`:

- DFT phonons: σ = 0.138;
- SevenNet phonons: σ = √(0.138² + s²), where s = 0.046 at 293 K and 0.040 at 77 K is the typical (68th
  percentile) difference between SevenNet and DFT for the same crystal;
- soft lattices (`soft_lattice`): at least 0.144.

The thickness is then known to within a factor exp(±σ): L₀.₈·exp(−σ) to L₀.₈·exp(+σ), about ±15 %. Rejections scale
roughly with thickness, so the `_lo` and `_hi` columns give the band of each rejection.

**Using it to compare two crystals.** Compare at the same temperature and on the same figure of merit. As a rule of
thumb, treat a difference as meaningful only if the two bands `[_lo, _hi]` do not overlap, which for two DFT
crystals means a difference in ln L₀.₈ of more than about 0.28 (thicknesses differing by more than about 30 %).
This rule treats the two errors as independent. Part of the model error is shared by similar crystals, so for two
crystals of similar chemistry the rule is conservative; for a probability rather than a yes or no, the pipeline
draws the errors by Monte Carlo (step 10, `rank_mc` in `configs/shortlist.yaml`, results in
`stage4_rank_probabilities.json`).

## Domain and practicality flags

Flags mark where the model was not validated and what makes a crystal impractical. They do not change any number.

| Flag | Set when | Effect |
|---|---|---|
| `soft_lattice` | mean phonon energy below 10 meV | outside the validated domain |
| `heavy_soft` | soft lattice with an element of Z ≥ 80 | outside the validated domain |
| `molecular` | molecular crystal (zero-dimensional bonding, no metal) or noble-gas solid | outside the validated domain |
| `untested_heavy` | element with Z ≥ 80 in a lattice that is not soft | kept in the comparison; heavier than any crystal in the validation |
| `mlip_unstable` | SevenNet phonons are dynamically unstable | kept; values are less reliable |
| `interp_only` | the DFT or SevenNet phonons are unstable only between commensurate q-points | kept; usually an interpolation artefact |
| `unscreenable` | the force-constant stability screen could not run | kept; stability not confirmed |
| `tierA_phase_mismatch` | the grown phase may not be the matched Materials Project entry | kept; check the phase |
| `not_off_bragg_single_crystal` | the material is not used as an off-Bragg single crystal | kept |
| `flag_Be`, `flag_Tl`, `flag_As`, `flag_Pb`, `flag_Bi` | contains that element (toxicity and handling) | practicality only |
| `hygroscopic` | belongs to a water-sensitive family (approximate rule) | practicality only |
| `melting_point_unknown`, `nucdata_incomplete` | the atlas holds no melting points; an element lacks nuclear data | informational |

`outside_domain` = `soft_lattice` or `heavy_soft` or `molecular`. `main_front_eligible` = a full model and not
`outside_domain`: the crystals of the main comparison.

## Validation against measured cross-sections

The model is checked against measured total cross-sections from EXFOR for 13 crystals: Al₂O₃, Si, SiO₂ (quartz), MgO,
Ge, Cu, Zn, BeO, graphite (C), Pb, Bi, Al and Be; a 14th, ZrO₂, is reported but not counted because the phase of its
samples is unknown. Crystals outside this list are not directly validated; the model uncertainty is what the
comparison on these crystals supports. For each dataset and model the statistic is the root-mean-square of ln(model/data) over 1–100 meV
(`r_d`), on logarithmic bins of 20 per decade. Single-crystal data are compared with the off-Bragg baseline at the
20th percentile of the data in each bin (the lower envelope beneath the Bragg spikes); polycrystal and powder data
with the full model including Bragg scattering, at the median. A crystal fails if its RMS exceeds 0.15 and also exceeds the
spread between its own datasets; the acceptance criterion, fixed before the computation, fails if 6 or more of the 13
crystals fail. The comparison was afterwards also restricted to single-crystal data, the only sample form the model
describes. The outcome and the per-crystal values are in `stage2_exfor.json` (key `crystals`) and
`stage2_kill2_D11.json` (the count on single-crystal data), with the per-dataset comparisons in
`stage2_exfor_benchmark.parquet` ([Results tables](results-tables.md)).

## Temperatures

Every crystal is computed at **293.15 K** (room temperature) and **77 K** (liquid-nitrogen cooling). Cooling
depopulates the phonons, so inelastic scattering and with it μ fall, and L₀.₈ grows; the rejections at the larger
thickness grow with it. The fast-neutron and photon cross sections are taken at 294 K for both temperatures (they
hardly depend on temperature); the moderator spectrum stays at 300 K. Other crystal temperatures require a rerun
([Extending](extending.md#3-change-the-temperature)); a series from 293 K to 77 K for the commercially grown shortlist
is in `stage4_temperature.json`.
