# Fab Access Simulator

A small, dependency-free companion to the Data-Rich, Insight-Poor? essay “What the Antibody Can Still Reach.” It asks when a longer-lived antibody can compensate for slower target capture, faster washout, or a smaller fraction of reachable epitopes. The simulation starts from free antibody in tissue, not serum.

## Run

Requires Python 3.10 or newer. From this repository’s directory:

```sh
python fab_access.py --scenario example.json --out results --sweep
python -m unittest discover -s tests -v
```

Open `results/comparison.svg` for the two-panel plot. `timecourses.csv` contains the plotted values, `summary.csv` gives the occupancy metrics, `inputs.json` records the exact inputs, and `sweep.csv` varies the variant’s association rate and accessible fraction independently. To model a matched concentration clamp, add `"constant_concentration_nM": 3` to the scenario JSON. The longer-exposure comparison then disappears by design.

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

Use nM for concentration and hours for time. `kon` is in nM⁻¹ h⁻¹; `koff` and `kloss` are in h⁻¹. The optional `bound_loss_per_h` represents a first-order loss of occupied sites with a fixed, replenished surface capacity; it is not an endosomal trafficking model. With constant exposure and `kloss = 0`, the occupancy ceiling is `f` and the equilibrium is `f × C / (C + koff/kon)` for positive `kon`. With varying exposure, the solver applies the exact constant-input update over each time step using the midpoint concentration.

This effective one-state model does not explicitly resolve two Fab arms, receptor diffusion or clustering, ligand shedding, rebinding, drug depletion, or target-mediated drug disposition. If cell-surface avidity is folded into an effective rate, that rate should be interpreted only under the antigen density and presentation used to estimate it. The model predicts occupancy, not neutralization, signaling, efficacy, or dose.

## Published context versus model assumptions

[Weise et al., bioRxiv, 2026, Fig. 8](https://www.biorxiv.org/content/10.64898/2026.11.02.763654v1.full) reported serum clearance values of 27.6 and 8.9 mL/day/kg for parental and inverted golimumab in IVIG-supplemented human-FcRn mice, and soluble TNF reporter IC₅₀ values of 38 and 28 pM. The preprint did not report paired cell-surface binding rates or free tissue concentration profiles for these formats. None of those reported values sets the simulator’s defaults; the example is synthetic. A threefold difference in serum clearance cannot be substituted for a threefold increase in tissue concentration or half-life.

[Shealy et al., mAbs, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3180089/) compared parental golimumab binding to soluble and transmembrane TNF using different assays. [Arnold et al., Nature Communications, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10199101/) used controlled antigen heights to show that membrane context can change effective antibody binding. These studies motivate the model’s distinction between exposure and access; they do not calibrate the variant.
