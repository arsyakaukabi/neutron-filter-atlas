# Glossary

What this page is for: one-line definitions of the technical terms used in this documentation, in alphabetical order.
Who should read it: anyone who meets an unfamiliar term on another page.

| Term | Definition |
|---|---|
| **1/E spectrum** | neutron flux per unit energy proportional to 1/E, the slowing-down spectrum of a moderator; the atlas uses it over 1 eV–100 keV |
| **Absorption** | capture of a neutron by a nucleus; for most elements its cross section grows linearly with wavelength (1/v law) |
| **Acceptance criterion** | the rule, fixed before the computation, that decides whether the model agrees with measured cross-sections: it fails if 6 or more of 13 crystals deviate |
| **Barn (b)** | unit of cross section, 10⁻²⁴ cm² |
| **Bragg scattering** | coherent elastic scattering by a crystal at wavelengths and orientations that satisfy Bragg's law; removed in the off-Bragg baseline |
| **Bragg dip** | a narrow drop of a single crystal's transmission where an orientation satisfies Bragg's law |
| **Commensurate q-point** | a phonon wavevector that fits the supercell exactly, where finite-displacement phonons are exact |
| **Commercially grown crystal** | a crystal offered for sale by a vendor, matched to a Materials Project entry; tier A |
| **Composition proxy** | thermal attenuation from free-element absorption and incoherent scattering only, without phonons |
| **Convex hull** | the lowest-energy combinations of compounds at each composition; distance above it (e_hull) measures thermodynamic stability |
| **Cross section (σ)** | probability of an interaction per atom, as an area (barn) |
| **Debye model** | a lattice-dynamics model with a single characteristic temperature θ_D for the phonon spectrum |
| **DFPT** | density-functional perturbation theory; one of the two Materials Project phonon workflows (`dfpt`) |
| **DFT** | density-functional theory, the electronic-structure method behind the Materials Project phonons |
| **DOS rule** | a crystal is dynamically unstable if more than 0.5 % of its phonon DOS lies below −0.5 THz |
| **Domain flag** | a flag that marks a crystal outside the validated range of the model (soft lattice, heavy soft lattice, molecular crystal) |
| **Dynamical stability** | absence of imaginary phonon frequencies; an unstable phonon calculation cannot be used |
| **Elastic Debye temperature** | θ_D computed from the elastic moduli (Anderson's formula) |
| **ENDF/B-VIII.0** | the evaluated nuclear-data library used for the fast-neutron and photon cross sections |
| **Epithermal neutrons** | neutrons above the thermal range, here 1 eV–100 keV |
| **EXFOR** | the international database of measured nuclear-reaction data, here total cross-sections |
| **Figure of merit** | a number that measures how good a filter is; in the atlas, the rejections at L₀.₈ |
| **Finite displacements** | computing phonons from the forces on atoms displaced by a small step in a supercell |
| **Force constants** | second derivatives of the energy with respect to atomic displacements, from which phonons follow |
| **Front** | see Pareto front |
| **Full phonon model** | thermal attenuation from per-element vibrational DOS of a phonon calculation (DFT or learned potential) |
| **Funnel** | the sequence of selections that decides which level of physics each crystal receives |
| **Fidelity** | how much physics a value includes (composition proxy, elastic Debye or full phonon model); also used in file and key names for the comparison of these levels |
| **Hot band** | wavelengths 0.3–1.0 Å, for hot-source instruments |
| **Incoherent scattering** | scattering without interference between atoms; elastic and inelastic parts |
| **Inelastic scattering** | scattering in which a neutron gains or loses energy by creating or absorbing phonons |
| **Kendall tau** | a rank correlation between two orderings of the same items (1 = identical order) |
| **L₀.₈** | the thickness at which the off-Bragg crystal transmits 80 % of the Maxwellian-weighted thermal band |
| **Learned interatomic potential** | a machine-learned model of the energy and forces of atoms, used here instead of DFT for phonons (SevenNet) |
| **Maxwellian** | the flux distribution of neutrons in thermal equilibrium with a moderator, here at 300 K |
| **MLIP** | machine-learned interatomic potential, written so in file, key and column names |
| **Model uncertainty (σ_model)** | the standard uncertainty of ln L₀.₈ from the comparison with measured cross-sections |
| **Mosaic spread** | the angular spread of crystal blocks in a real single crystal |
| **NCrystal** | the library that computes thermal-neutron cross sections of crystals from their structure and dynamics |
| **NCMAT** | NCrystal's material-file format |
| **Nieman ratio** | epithermal attenuation at 1 eV divided by thermal attenuation at 1.8 Å |
| **Off-Bragg** | describing a single crystal oriented so that no Bragg reflection is excited in the band of interest |
| **OpenMC** | the Monte Carlo transport code used for slab transport, activation and the nuclear-data extraction |
| **Pareto front** | the crystals not beaten on both of two figures of merit by any other crystal |
| **pheasy** | the finite-displacement phonon workflow of the Materials Project (`pheasy` in file and column names) |
| **Phonon** | a quantum of lattice vibration |
| **Phonon DOS** | density of states of the phonons as a function of frequency or energy |
| **Practicality flag** | a flag about toxicity, water sensitivity or handling; it changes no number |
| **Proxy bound** | the selection that removes crystals whose composition-only transmission through 1 cm is below 0.5 |
| **Rejection** | −ln of the transmission of an unwanted component (epithermal or fission neutrons, photons) |
| **SevenNet** | the learned interatomic potential used for the atlas phonons (model SevenNet-MF-ompa) |
| **Shortlist** | the crystals examined further (transport, orientation, temperature, practicality) |
| **SIMEX** | simulation-extrapolation, a method that estimates the effect of noise on a statistic; used in a reproducibility check |
| **Soft lattice** | a crystal whose mean phonon energy is below 10 meV |
| **Supercell** | a repetition of the unit cell used for finite-displacement phonons |
| **`stage0_` to `stage4_`** | prefixes of the result file names, grouping the files by the part of the computation that produced them; see [Data](data.md#what-the-results-archive-holds) |
| **Thermal band** | wavelengths 0.5–3.0 Å, over which the thermal transmission is averaged |
| **Thermal transmission (T_th)** | the Maxwellian-weighted average transmission over the thermal band |
| **Tier A, tier B** | commercially grown crystals (A) and all other candidates (B) |
| **Triviality test** | the test whether the composition proxy orders crystals as well as the full model |
| **Uncollided** | crossing a slab without any interaction; the atlas rejections count only uncollided particles |
| **VDOS** | vibrational density of states per element, the input NCrystal needs for inelastic scattering |
| **Watt spectrum** | the energy spectrum of fission neutrons, S(E) ∝ exp(−E/a) sinh(√(bE)) |
