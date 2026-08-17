#!/usr/bin/env python3
"""Regenerate the pinned Tensorpac formula fixture for WS9-16."""

import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import tensorpac
from tensorpac.methods import mean_vector_length, modulation_index


HERE = Path(__file__).resolve().parent
OUT = HERE / "fixtures" / "tensorpac-0.6.5"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


n = 512
index = np.arange(n, dtype=float)
phase = -math.pi + 2 * math.pi * (index + 0.37) / n
amplitude = (
    1.1
    + 0.4 * np.cos(phase - 0.31)
    + 0.15 * np.sin(3 * phase)
    + 0.01 * (index % 7)
)

with (OUT / "input.csv").open("w", newline="", encoding="ascii") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["phase", "amplitude"])
    writer.writerows(zip(phase, amplitude))

tort = float(modulation_index(phase[None, :], amplitude[None, :], 18).item())
canolty = float(mean_vector_length(phase[None, :], amplitude[None, :]).item())
ozkurt = float(
    abs(np.sum(amplitude * np.exp(1j * phase)))
    / math.sqrt(n * np.sum(amplitude**2))
)

with (OUT / "expected.csv").open("w", newline="", encoding="ascii") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["metric", "value", "reference"])
    writer.writerow(["tort", format(tort, ".17g"), "tensorpac.methods.modulation_index"])
    writer.writerow(["canolty", format(canolty, ".17g"),
                     "tensorpac.methods.mean_vector_length"])
    writer.writerow(["ozkurt", format(ozkurt, ".17g"),
                     "published_formula_independent_numpy"])

generator_hash = sha256(Path(__file__).resolve())
input_hash = sha256(OUT / "input.csv")
expected_hash = sha256(OUT / "expected.csv")
manifest = "\n".join([
    "Fixture-ID: ws9-16-tensorpac-formulas",
    "WP-ID: WS9-16",
    "Reference-Tool: Tensorpac",
    f"Reference-Version: {tensorpac.__version__}",
    "Reference-URL: https://etiennecmb.github.io/tensorpac/",
    "Generator: generate-ws9-16-tensorpac.py",
    f"Generator-SHA256: {generator_hash}",
    f"Input-SHA256: {input_hash}",
    f"Expected-SHA256: {expected_hash}",
    "Runtime: CPython 3.14; NumPy 2.5.1",
    "Seed: deterministic closed-form input; no RNG",
    "Direction-Convention: phase first, amplitude second; upper-tail PAC",
    "Generated-UTC: 2026-07-29T00:00:00Z",
    ("License-Note: numeric outputs only; Tensorpac source is not copied; "
     "Tensorpac is BSD-3-Clause"),
    ("Ozkurt-Note: Tensorpac norm_direct_pac uses a different normalization; "
     "the Ozkurt fixture uses the published package formula in independent NumPy"),
    "",
])
(OUT / "manifest.dcf").write_text(manifest, encoding="ascii")

covered = [
    OUT / "input.csv",
    OUT / "expected.csv",
    OUT / "manifest.dcf",
    Path(__file__).resolve(),
]
lines = []
for path in covered:
    display = path.name if path.parent == OUT else f"../../{path.name}"
    lines.append(f"{sha256(path)}  {display}")
(OUT / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="ascii")

print(f"wrote {OUT}")
