#!/usr/bin/env python3
"""
Calculate total and Burgers-vector-resolved dislocation densities from an
OVITO/LAMMPS trajectory sequence using the Dislocation Analysis (DXA) modifier.

Example
-------
python dislocation_density.py ^
  --input "F:/TZ/Ni-Au/3/Shuju/Indentation.*" ^
  --output "dislocation_density.tsv" ^
  --frame0-step 51000 ^
  --step-interval 1000 ^
  --min-step 100000 ^
  --max-step 151000 ^
  --density-unit m-2
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
from ovito.io import import_file
from ovito.modifiers import DislocationAnalysisModifier


DISLOCATION_TYPES = {
    "Other": None,
    "Perfect": "1/2<110>",
    "Shockley": "1/6<112>",
    "Stair-rod": "1/6<110>",
    "Hirth": "1/3<100>",
    "Frank": "1/3<111>",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate dislocation densities from an OVITO trajectory."
    )
    parser.add_argument(
        "--input",
        required=True,
        help='Input trajectory pattern, for example "Indentation.*".',
    )
    parser.add_argument(
        "--output",
        default="dislocation_density.tsv",
        help="Output TSV file. Existing files are overwritten.",
    )
    parser.add_argument(
        "--frame0-step",
        type=int,
        default=51000,
        help="LAMMPS step corresponding to frame 0.",
    )
    parser.add_argument(
        "--step-interval",
        type=int,
        default=1000,
        help="LAMMPS step interval between consecutive frames.",
    )
    parser.add_argument(
        "--min-step",
        type=int,
        default=100000,
        help="First LAMMPS step to include.",
    )
    parser.add_argument(
        "--max-step",
        type=int,
        default=151000,
        help="Last LAMMPS step to include.",
    )
    parser.add_argument(
        "--density-unit",
        choices=("A-2", "m-2"),
        default="A-2",
        help="Output density unit. OVITO length/volume gives A^-2 when coordinates are in angstroms.",
    )
    return parser.parse_args()


def get_cell_volume(data) -> float:
    """Return the simulation-cell volume in the trajectory's length unit cubed."""
    volume = data.attributes.get("DislocationAnalysis.cell_volume")
    if volume is not None:
        return float(volume)

    if data.cell is None:
        raise RuntimeError("The frame does not contain simulation-cell information.")

    # OVITO's cell matrix stores the three cell vectors in the first 3 columns.
    matrix = np.asarray(data.cell.matrix)
    if matrix.shape[0] < 3 or matrix.shape[1] < 3:
        raise RuntimeError(f"Unexpected simulation-cell matrix shape: {matrix.shape}")

    a, b, c = matrix[:, 0], matrix[:, 1], matrix[:, 2]
    volume = abs(float(np.dot(a, np.cross(b, c))))
    if volume <= 0.0:
        raise RuntimeError(f"Invalid simulation-cell volume: {volume}")
    return volume


def find_dislocation_length(attributes, label: str, burgers: str | None) -> float:
    """Find the line length for one dislocation type in OVITO attributes."""
    # Common OVITO attribute naming conventions.
    candidates = []
    if burgers is not None:
        candidates.extend(
            [
                f"DislocationAnalysis.length.{burgers}",
                f"DislocationAnalysis.length.{label}",
            ]
        )
    else:
        candidates.extend(
            [
                "DislocationAnalysis.length.other",
                "DislocationAnalysis.length.Other",
            ]
        )

    for key in candidates:
        if key in attributes:
            return float(attributes[key])

    # Fallback for small naming differences between OVITO versions.
    label_lower = label.lower()
    burgers_lower = burgers.lower() if burgers else None
    for key, value in attributes.items():
        key_lower = key.lower()
        if "dislocationanalysis" not in key_lower or "length" not in key_lower:
            continue
        if label_lower in key_lower or (burgers_lower and burgers_lower in key_lower):
            return float(value)

    return 0.0


def convert_density(value_a2: float, unit: str) -> float:
    """Convert A^-2 to the requested output unit."""
    if unit == "m-2":
        return value_a2 * 1.0e20
    return value_a2


def main() -> None:
    args = parse_args()

    pipeline = import_file(args.input)
    dxa = DislocationAnalysisModifier()
    dxa.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
    pipeline.modifiers.append(dxa)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    density_label = "m^-2" if args.density_unit == "m-2" else "A^-2"
    column_names = ["Step", f"Total ({density_label})"]
    column_names.extend(f"{name} ({density_label})" for name in DISLOCATION_TYPES)

    processed_frames = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(column_names) + "\n")

        for frame in range(pipeline.source.num_frames):
            current_step = args.frame0_step + frame * args.step_interval

            if current_step < args.min_step:
                continue
            if current_step > args.max_step:
                break

            data = pipeline.compute(frame)
            cell_volume = get_cell_volume(data)

            total_length = float(
                data.attributes.get("DislocationAnalysis.total_line_length", 0.0)
            )
            densities: Dict[str, float] = {
                "Total": convert_density(total_length / cell_volume, args.density_unit)
            }

            for label, burgers in DISLOCATION_TYPES.items():
                line_length = find_dislocation_length(
                    data.attributes, label=label, burgers=burgers
                )
                densities[label] = convert_density(
                    line_length / cell_volume, args.density_unit
                )

            values = [str(current_step), f"{densities['Total']:.10e}"]
            values.extend(f"{densities[name]:.10e}" for name in DISLOCATION_TYPES)
            handle.write("\t".join(values) + "\n")
            handle.flush()

            processed_frames += 1
            print(f"Step {current_step} processed.")

    if processed_frames == 0:
        raise RuntimeError(
            "No frames were written. Check --frame0-step, --step-interval, "
            "--min-step, and --max-step."
        )

    print(f"Completed. Results saved to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
