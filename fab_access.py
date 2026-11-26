"""Toy model for separating tissue exposure, Fab access, and occupancy."""
import argparse
from bisect import bisect_right
import csv
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path
from xml.sax.saxutils import escape


@dataclass(frozen=True)
class Exposure:
    initial_nM: float
    half_life_h: float

    def __post_init__(self):
        if not all(map(math.isfinite, (self.initial_nM, self.half_life_h))):
            raise ValueError("Exposure parameters must be finite.")
        if self.initial_nM < 0 or self.half_life_h <= 0:
            raise ValueError("Initial concentration must be nonnegative; half-life must be positive.")

    def at(self, hours):
        return self.initial_nM * 2.0 ** (-hours / self.half_life_h)

    @property
    def knots(self):
        return ()


@dataclass(frozen=True)
class Profile:
    points: tuple

    def __post_init__(self):
        if not 2 <= len(self.points) <= 2000 or any(len(pair) != 2 for pair in self.points):
            raise ValueError("Provide 2 to 2,000 time/concentration pairs.")
        if any(not math.isfinite(x) or x < 0 for pair in self.points for x in pair):
            raise ValueError("Profile times and concentrations must be finite and nonnegative.")
        times = tuple(t for t, _ in self.points)
        if times[0] != 0 or any(y <= x for x, y in zip(times, times[1:])):
            raise ValueError("Profile times must start at zero and strictly increase.")
        object.__setattr__(self, "_knots", times)

    @property
    def knots(self):
        return self._knots

    def at(self, hours):
        i = min(len(self.points) - 2, bisect_right(self.knots, hours) - 1)
        (a, y), (b, z) = self.points[i:i+2]
        return y + (z - y) * ((min(b, hours) - a) / (b - a))


def parse_exposure(p):
    if "points" in p:
        return Profile(tuple(tuple(row) for row in p["points"]))
    return Exposure(**p)


@dataclass(frozen=True)
class Binding:
    kon_per_nM_h: float
    koff_per_h: float
    accessible_fraction: float = 1.0
    bound_loss_per_h: float = 0.0

    def __post_init__(self):
        rates = (self.kon_per_nM_h, self.koff_per_h, self.bound_loss_per_h)
        if any(not math.isfinite(x) or x < 0 for x in rates):
            raise ValueError("Binding rates must be finite and nonnegative.")
        if not 0 <= self.accessible_fraction <= 1:
            raise ValueError("Accessible fraction must be in [0, 1].")


def step(bound, concentration, binding, dt):
    """Exact occupancy step at a constant concentration."""
    association = binding.kon_per_nM_h * concentration
    total = association + binding.koff_per_h + binding.bound_loss_per_h
    if total == 0:
        return bound
    steady = binding.accessible_fraction * (association / total)
    return bound + (steady - bound) * -math.expm1(-total * dt)


def simulate(exposure, binding, duration_h, dt_h, constant_nM=None, extra_knots=()):
    if not all(map(math.isfinite, (duration_h, dt_h))) or min(duration_h, dt_h) <= 0:
        raise ValueError("Duration and step must be finite and positive.")
    if constant_nM is not None and (not math.isfinite(constant_nM) or constant_nM < 0):
        raise ValueError("Constant concentration must be finite and nonnegative.")
    ratio = duration_h / dt_h
    if ratio > 1_000_000:
        raise ValueError("Step size creates more than one million intervals.")
    count = max(1, math.ceil(ratio))
    dt = duration_h / count
    if exposure.knots and exposure.knots[-1] < duration_h:
        raise ValueError("Profile must cover the full simulated duration.")
    grid = {i * dt for i in range(count)} | {duration_h}
    if constant_nM is None:
        grid.update(t for t in (*exposure.knots, *extra_knots) if 0 < t < duration_h)
    concentration = exposure.at if constant_nM is None else lambda _: constant_nM
    rows = [(0.0, concentration(0), 0.0)]
    for time in sorted(grid)[1:]:
        prior, _, bound = rows[-1]
        bound = step(bound, concentration((prior + time) / 2), binding, time - prior)
        rows.append((time, concentration(time), bound))
    return rows


def summarize(rows, threshold):
    if not 0 < threshold <= 1:
        raise ValueError("Occupancy threshold must be in (0, 1].")
    time_above = exposure_area = occupancy_area = 0.0
    for (t0, c0, b0), (t1, c1, b1) in zip(rows, rows[1:]):
        dt = t1 - t0
        exposure_area += (c0 + c1) / 2 * dt
        occupancy_area += (b0 + b1) / 2 * dt
        if min(b0, b1) >= threshold:
            time_above += dt
        elif max(b0, b1) > threshold:
            time_above += dt * (max(b0, b1) - threshold) / abs(b1 - b0)
    return {
        "tissue_exposure_auc_nM_h": exposure_area,
        "occupied_fraction_auc_h": occupancy_area,
        "mean_occupied_fraction": occupancy_area / rows[-1][0],
        "peak_occupied_fraction": max(b for _, _, b in rows),
        "end_occupied_fraction": rows[-1][2],
        "hours_at_or_above_threshold": time_above,
    }


def compare(cfg):
    reference_exp = parse_exposure(cfg["reference_exposure"])
    variant_exp = parse_exposure(cfg["variant_exposure"])
    reference_binding = Binding(**cfg["reference_binding"])
    variant_binding = Binding(**cfg["variant_binding"])
    arguments = (
        cfg["duration_h"], cfg["dt_h"], cfg.get("constant_concentration_nM"),
        (*reference_exp.knots, *variant_exp.knots),
    )
    return {
        "reference": simulate(reference_exp, reference_binding, *arguments),
        "variant_same_exposure": simulate(reference_exp, variant_binding, *arguments),
        "exposure_only": simulate(variant_exp, reference_binding, *arguments),
        "variant_own_exposure": simulate(variant_exp, variant_binding, *arguments),
    }


def sweep(cfg):
    """Hold variant exposure fixed and vary access and association independently."""
    exp = parse_exposure(cfg["variant_exposure"])
    ref = parse_exposure(cfg["reference_exposure"])
    binding = Binding(**cfg["variant_binding"])
    metrics = compare(cfg)
    reference = summarize(metrics["reference"], cfg["occupancy_threshold"])
    for cap in (0.25, 0.50, 0.65, 0.80, 0.90, 1.00):
        for multiplier in (0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00):
            candidate = replace(
                binding, accessible_fraction=cap,
                kon_per_nM_h=binding.kon_per_nM_h * multiplier,
            )
            rows = simulate(exp, candidate, cfg["duration_h"], cfg["dt_h"],
                            cfg.get("constant_concentration_nM"), ref.knots)
            summary = summarize(rows, cfg["occupancy_threshold"])
            yield {
                "accessible_fraction": cap,
                "kon_multiplier_vs_variant": multiplier,
                "mean_occupied_fraction": summary["mean_occupied_fraction"],
                "delta_mean_occupancy_vs_reference": summary["mean_occupied_fraction"] - reference["mean_occupied_fraction"],
                "hours_at_or_above_threshold": summary["hours_at_or_above_threshold"],
                "delta_hours_vs_reference": summary["hours_at_or_above_threshold"] - reference["hours_at_or_above_threshold"],
            }


def svg_chart(series, threshold, scenario):
    colors = ("#123955", "#ab501d", "#667584", "#007b69")
    dashes = ("", "8 4", "3 5", "")
    names = ("Reference", "Variant, same exposure", "Exposure alone", "Variant, own exposure")
    duration = next(iter(series.values()))[-1][0]
    concentration_max = max(1e-9, max(c for rows in series.values() for _, c, _ in rows))
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="590" viewBox="0 0 1000 590" role="img" aria-labelledby="title desc">',
           '<title id="title">Fab access and exposure simulation</title>',
           '<desc id="desc">Two charts compare free tissue exposure and target occupancy in four hypothetical scenarios.</desc>',
           '<rect width="1000" height="590" fill="#fcfbf7"/>',
           '<g font-family="system-ui, sans-serif" fill="#183242">',
           '<text x="65" y="38" font-size="25" font-weight="700">Fab access and exposure</text>',
           f'<text x="65" y="63" font-size="14">{escape(scenario)} · Illustrative parameters, not fitted data</text>']
    for j, name in enumerate(names):
        x = 65 + j * 226
        out += [f'<line x1="{x}" x2="{x+25}" y1="542" y2="542" stroke="{colors[j]}" stroke-width="3" stroke-dasharray="{dashes[j]}"/>',
                f'<text x="{x+33}" y="547" font-size="13">{name}</text>']
    for j, (label, maximum) in enumerate((("Free antibody in tissue (nM)", concentration_max), ("Occupied surface target (%)", 1.0))):
        left, width, top, height = 70 + 485 * j, 380, 145, 315
        out.append(f'<text x="{left}" y="{top-26}" font-size="17" font-weight="650">{label}</text>')
        for i in range(5):
            y = top + height - i * height / 4
            value = maximum * i / 4 * (100 if j else 1)
            out += [f'<line x1="{left}" x2="{left+width}" y1="{y}" y2="{y}" stroke="#dce3e3"/>',
                    f'<text x="{left-8}" y="{y+4}" text-anchor="end" font-size="12">{value:.3g}</text>']
            x = left + width * i / 4
            out.append(f'<text x="{x}" y="{top+height+22}" text-anchor="middle" font-size="12">{duration*i/4:.3g}</text>')
        out.append(f'<text x="{left+width/2}" y="{top+height+45}" text-anchor="middle" font-size="13">Time (hours)</text>')
        if j:
            y = top + height * (1-threshold)
            out.append(f'<line x1="{left}" x2="{left+width}" y1="{y}" y2="{y}" stroke="#8b9194" stroke-dasharray="2 5"/>')
        else:
            out.append(f'<text x="{left+width}" y="{top-8}" text-anchor="end" font-size="11">Pairs overlap by design</text>')
        for i, rows in enumerate(series.values()):
            if not j and i in (1, 2):  # Equal exposures would overdraw their matching lines.
                continue
            stride = max(1, math.ceil(len(rows) / 1400))
            samples = rows[::stride]
            if samples[-1] != rows[-1]:
                samples.append(rows[-1])
            points = " ".join(f"{left + width*t/duration:.2f},{top+height*(1-v/maximum):.2f}"
                              for t, c, b in samples for v in [b if j else c])
            out.append(f'<polyline points="{points}" fill="none" stroke="{colors[i]}" stroke-width="2.7" stroke-dasharray="{dashes[i]}"/>')
    out += ['</g>', '</svg>']
    return "\n".join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scenario", type=Path, default=Path(__file__).with_name("example.json"))
    p.add_argument("--out", type=Path, default=Path("results"))
    p.add_argument("--sweep", action="store_true", help="Write a grid of access and association scenarios.")
    args = p.parse_args()
    try:
        cfg = json.loads(args.scenario.read_text())
        series = compare(cfg)
        threshold = cfg["occupancy_threshold"]
        summaries = [{"condition": key, **summarize(rows, threshold)} for key, rows in series.items()]
        if args.sweep and cfg["duration_h"] / cfg["dt_h"] > 50_000:
            raise ValueError("Use 50,000 or fewer intervals per run for the sweep.")
        args.out.mkdir(parents=True, exist_ok=True)
        with (args.out / "timecourses.csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(("condition", "time_h", "free_tissue_nM", "fraction_total_target_occupied"))
            for key, rows in series.items():
                writer.writerows((key, *row) for row in rows)
        blocks = [("summary.csv", summaries)]
        if args.sweep:
            blocks.append(("sweep.csv", list(sweep(cfg))))
        for filename, data in blocks:
            with (args.out / filename).open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(data[0]))
                writer.writeheader()
                writer.writerows(data)
        (args.out / "inputs.json").write_text(json.dumps(cfg, indent=2) + "\n")
        (args.out / "comparison.svg").write_text(svg_chart(series, threshold, cfg.get("label", "Custom scenario")))
        print(f"{'Condition':27}  {'Mean occupancy':>14}  {'Hours at/above threshold':>24}")
        for row in summaries:
            print(f"{row['condition']:27}  {row['mean_occupied_fraction']:14.1%}  {row['hours_at_or_above_threshold']:24.2f}")
        print(f"\nWrote CSV, input snapshot, and SVG files to {args.out.resolve()}")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        p.error(str(exc))


if __name__ == "__main__":
    main()
