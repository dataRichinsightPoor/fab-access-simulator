# Fab Access Simulator

A small, dependency-free companion to the Data-Rich, Insight-Poor? essay “What the Antibody Can Still Reach.” It asks when a longer-lived antibody can compensate for slower target capture, faster washout, or a smaller fraction of reachable epitopes. The simulation starts from free antibody in tissue, not serum.

## Run in your browser

Open the [Fab Access Simulator](https://datarichinsightpoor.github.io/fab-access-simulator/). Begin with a preset or change both formats independently. Give each format an exponential tissue-exposure profile or paste a custom `time_h,free_tissue_nM` CSV with increasing times starting at zero and ending at or after the simulation duration. For a pulse and washout, place two closely spaced points on either side of the transition. The browser interpolates between input points, so this is a short ramp rather than a mathematical discontinuity. All calculation and file handling happen locally; the app sends no scenario data to a backend.

The four curves separate exposure from binding; the sensitivity grid compares time at or above the user-set threshold across accessible fractions and variant association rates. Export a scenario as JSON for reuse in the browser or Python, and export the computed trajectories or the sensitivity grid as CSV.

## Run locally

Requires Python 3.10 or newer. From this repository’s directory:

```sh
python fab_access.py --scenario example.json --out results --sweep
python -m unittest discover -s tests -v
node --test tests/*.test.mjs
```

Open `results/comparison.svg` for the two-panel plot. `timecourses.csv` contains the plotted values, `summary.csv` gives the occupancy metrics, `inputs.json` records the exact inputs, and `sweep.csv` varies the variant’s association rate and accessible fraction independently. To model a matched concentration clamp, give both formats the same constant custom profile or add `"constant_concentration_nM": 3` to the scenario JSON. The longer-exposure comparison then disappears by design. To run the browser locally, `python -m http.server 8000 --directory docs` serves the app on `http://localhost:8000`.

The [browser solver](docs/model.mjs) and [Python solver](fab_access.py) are separate implementations. The [cross-runtime tests](tests/model.test.mjs) compare the default scenario sample-by-sample, check analytic constant-exposure and washout cases, and test exposure knots, hard access limits, and threshold integration. [Python tests](tests/test_model.py) also check numerical convergence and shared time grids. Node.js 20+ is recommended for running JS tests; it is not needed to use the browser interface.

## What the four curves distinguish

| Condition | Tissue exposure | Binding |
|---|---|---|
| Reference | Reference profile | Reference parameters |
| Variant, same exposure | Reference profile | Variant parameters |
| Exposure alone | Variant profile | Reference parameters |
| Variant, own exposure | Variant profile | Variant parameters |

The default scenario deliberately includes a trade-off. The variant has a longer assumed tissue half-life but weaker association and a lower access ceiling. Change one factor at a time before trying to explain the combined curves. “Hours at or above threshold” uses linear interpolation between time samples; the threshold is a modeling choice, not a biological efficacy cutoff.

Across the default 72-hour window, the variant at its own exposure averages 44.2% occupancy versus 42.6% for the reference, yet spends only 27.2 versus 29.5 hours at or above 50% occupancy. The opposite ranking is intentional: integrating occupancy and counting time above a threshold are different questions.

## Model

The model imposes a monoexponential free tissue concentration `C(t) = C0 × 2^(-t / t½)`, or a constant clamp. It tracks `B`, the fraction of all surface sites occupied, with an accessible fraction `f`:

```text
dB/dt = kon × C(t) × (f − B) − (koff + kloss) × B
B(0)  = 0
```

Use nM for concentration and hours for time. `kon` is in nM⁻¹ h⁻¹; `koff` and `kloss` are in h⁻¹. The optional `bound_loss_per_h` represents a first-order loss of occupied sites with a fixed, replenished surface capacity; it is not an endosomal trafficking model. For a step of length `h`, set `a = kon × C(midpoint)` and `r = a + koff + kloss`. When `r > 0`, the update is `B_next = B_inf + (B − B_inf) × exp(−r × h)`, with `B_inf = f × a/r`; when `r = 0`, occupancy stays unchanged. Profile knot times are included on a common time grid so comparisons respect changes in slope. Mean occupancy is the trapezoidal integral divided by the observation window. With varying concentration, rerun with a smaller maximum step and check that the metrics converge.

This effective one-state model does not explicitly resolve two Fab arms, receptor diffusion or clustering, ligand shedding, rebinding, drug depletion, or target-mediated drug disposition. If cell-surface avidity is folded into an effective rate, that rate should be interpreted only under the antigen density and presentation used to estimate it. The model predicts occupancy, not neutralization, signaling, efficacy, or dose.

## Published context versus model assumptions

[Weise et al., bioRxiv, 2026, Fig. 8](https://www.biorxiv.org/content/10.64898/2026.11.02.763654v1.full) reported serum clearance values of 27.6 and 8.9 mL/day/kg for parental and inverted golimumab in IVIG-supplemented human-FcRn mice, and soluble TNF reporter IC₅₀ values of 38 and 28 pM. The preprint did not report paired cell-surface binding rates or free tissue concentration profiles for these formats. None of those reported values sets the simulator’s defaults; the example is synthetic. A threefold difference in serum clearance cannot be substituted for a threefold increase in tissue concentration or half-life.

[Shealy et al., mAbs, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3180089/) compared parental golimumab binding to soluble and transmembrane TNF using different assays. [Arnold et al., Nature Communications, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10199101/) used controlled antigen heights to show that membrane context can change effective antibody binding. These studies motivate the model’s distinction between exposure and access; they do not calibrate the variant.

## Structural accuracy and limits

### Inspect the coordinate-based geometry

The geometry panel projects one deposited TNF trimer (chains A/B/C) and the Fv contacting chain A (light F, heavy I) from [PDB 5YOY](https://www.rcsb.org/structure/5YOY). Gold markers use the epitope residues reported by [Ono et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5980524/) in mature TNF numbering. The original structure stays rigid as it tilts relative to the hypothetical membrane. Change the origin height, tilt, and size and placement of an optional spherical domain envelope. The reference frame, pose and plane-overlap calculation are in [`geometry.mjs`](docs/geometry.mjs) and [`scripts/build-5yoy.mjs`](scripts/build-5yoy.mjs).

For each deposited heavy atom with coordinates `(x,y,z)` in the display frame, the signed distance to the hypothetical plane at `y=0` is `height + x*sin(tilt) + y*cos(tilt) − 0.17 nm`; the minimum determines the padded clearance of each component. For the sphere, clearance is its center height minus its radius. The minimum origin height is the lift needed from zero to give all displayed components nonnegative clearance. The sphere only tests the plane, not collisions with the TNF or Fab. It is not a molecular model of Fc.

The bundled subset is derived from public archival coordinates, separately credited to the PDB depositors. To reproduce it, download the [5YOY PDB file](https://files.rcsb.org/download/5YOY.pdb) and run `node scripts/build-5yoy.mjs 5YOY.pdb > docs/5yoy-data.mjs`. The generated module records the input’s SHA-256.

### Do not fit occupancy to the geometry

The kinetic model has no molecular coordinates, docking scores, explicit Fc or Fab linkers, glycocalyx, or computed steric footprint. The structural panel does not set `f` or any rate. Treat `f` as an effective accessible fraction over the experimental window, not a structural count of exposed sites. A lower apparent association rate and a lower `f` can both suppress the early signal; for an initially empty surface their contribution to initial slope is a product. A single trace generally cannot separate them. Matched concentration series, an independently assessed surface target level, and controls for target height and ligand presentation are needed to test those interpretations.

[PDB 5YOY](https://www.rcsb.org/structure/5YOY) resolves the golimumab Fv with soluble TNF. It provides an epitope footprint but no coordinates for the intact inverted antibody, its C-terminal Fab connections, or a cell membrane. A clear pose in the panel therefore does not mean it is populated or accessible in cells, and an overlap at an assumed tilt does not establish a structural penalty for golimumab.

## License

[MIT](LICENSE). The license covers this repository’s software and accompanying documentation, not the cited publications or deposited structures.
