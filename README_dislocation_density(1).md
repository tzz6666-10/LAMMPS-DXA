# OVITO DXA dislocation-density post-processing script

## Purpose

`dislocation_density.py` calculates the total dislocation density and the
densities of the following FCC dislocation types from an OVITO-readable
LAMMPS trajectory sequence:

- Perfect: `1/2<110>`
- Shockley: `1/6<112>`
- Stair-rod: `1/6<110>`
- Hirth: `1/3<100>`
- Frank: `1/3<111>`
- Other

The script applies OVITO's **Dislocation Analysis (DXA)** modifier with the
input crystal structure explicitly set to FCC.

## Requirements

- Python 3.9 or later
- OVITO Python module
- NumPy

Install the Python packages according to the OVITO documentation for your
operating system and Python environment.

## Default folders

- Input trajectory: `F:\\TZ\\Ni-Au\\3\\Shuju\\Indentation.*`
- Output file: `E:\\NMYH-denstiy-3\\dislocation_density_final.txt`

These defaults preserve the original folders used in the calculation. They can still be overridden with `--input` and `--output`.

## Example

```bash
python dislocation_density.py \
  --input "F:\\TZ\\Ni-Au\\3\\Shuju\\Indentation.*" \
  --output "E:\\NMYH-denstiy-3\\dislocation_density_final.txt" \
  --frame0-step 51000 \
  --step-interval 1000 \
  --min-step 100000 \
  --max-step 151000 \
  --density-unit A-2
```

On Windows Command Prompt, use `^` instead of `\` for line continuation, or
write the command on one line.

## Arguments

- `--input`: OVITO-readable file sequence or wildcard pattern. The default path is `F:\\TZ\\Ni-Au\\3\\Shuju\\Indentation.*`.
- `--output`: Output TSV path. The default path is `E:\\NMYH-denstiy-3\\dislocation_density_final.txt`. Existing files are overwritten.
- `--frame0-step`: LAMMPS step corresponding to the first imported frame.
- `--step-interval`: Step interval between consecutive frames.
- `--min-step`: First step included in the output.
- `--max-step`: Last step included in the output.
- `--density-unit`: `A-2` or `m-2`.

## Units

DXA reports line length divided by cell volume. If the trajectory coordinates
are in angstroms, the raw density unit is `Å^-2`. The script can convert this
to `m^-2` using:

`1 Å^-2 = 1e20 m^-2`.

## Output

The output is a tab-separated text file containing the LAMMPS step, total
dislocation density, and Burgers-vector-resolved dislocation densities.

## Notes

- Confirm that `--frame0-step` and `--step-interval` match the dump settings
  used in the LAMMPS simulation.
- The script overwrites the output file to avoid duplicated rows when it is
  rerun.
- Attribute names may vary slightly between OVITO versions. The script checks
  common names and includes a fallback keyword match.
