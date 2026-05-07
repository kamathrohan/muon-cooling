# muon-cooling

Python toolkit for designing muon ionization cooling channels and rendering them to BDSIM `.gmad` files.

---

## Contents

1. [Overview](#overview)
2. [Directory structure](#directory-structure)
3. [Quickstart](#quickstart)
4. [Physics](#physics)
   - [Element classes](#element-classes)
   - [Beamline and MuonCoolingChannel](#beamline-and-muoncoolingchannel)
   - [Builder functions](#builder-functions)
   - [DataFrame introspection](#dataframe-introspection)
5. [Rendering](#rendering)
   - [render_gmad](#render_gmad)
   - [generate_sampler_lines](#generate_sampler_lines)
   - [generate_beam_block](#generate_beam_block)
6. [Beam file translation](#beam-file-translation)
7. [JSON-driven pipeline](#json-driven-pipeline)

---

## Overview

A muon ionization cooling channel consists of repeating cells containing solenoid coils (focusing), wedge/cylinder absorbers (cooling), RF cavities (re-acceleration), and dipoles (dispersion). This library models those elements in Python, builds periodic channels programmatically, and renders them to BDSIM `.gmad` via a Jinja2 template.

**Main workflow:**

```
MuonCoolingChannel → build_*_beamline() → render_gmad() → .gmad
```

---

## Directory structure

```
muon-cooling/
├── src/
│   ├── physics/          # element classes, beamline, builders, translation
│   ├── rendering/        # GMAD rendering and Jinja2 templating
│   ├── pipeline.py       # config-driven build (JSON → gmad)
│   └── run.py            # config-driven build + simulate + analyse
├── analyse/              # post-processing scripts and notebooks
├── config/               # Jinja2 template and JSON config
├── scripts/              # HPC job scripts
├── output/               # simulation artefacts (gitignored)
└── quickstart.py         # working example
```

---

## Quickstart

```python
from src import (
    MuonCoolingChannel,
    build_coil_beamline, build_dipole_beamline,
    build_absorber_beamline, build_rf_beamline,
    render_gmad,
)

channel = MuonCoolingChannel(
    n_cells=10, cell_length=1.0, total_length=170.0,
    total_width=0.8, reference_momentum=200.0, on_axis_tolerance=2e-5,
)

coil_1 = dict(name='Coil1', r_in=0.285, r_thick=0.070,
              N_pancakes=17, L_pancake=0.012, L_spacing=0.004,
              currDensity=328.43, N_sheets=5)
coil_2 = dict(name='Coil2', r_in=0.185, r_thick=0.060,
              N_pancakes=5, L_pancake=0.012, L_spacing=0.003,
              currDensity=300.0, N_sheets=5)
dipole  = dict(name='Dipole', field_strength=0.2,
               aperture=0.2, length_z=0.1, enge_coefficient=5.5)
absorber = dict(name='Absorber', absorber_type='wedge',
                material='G4_LITHIUM_HYDRIDE',
                wedge_opening_angle=0.1745, wedge_height=0.24,
                wedge_apex_to_base=0.1134256)
rf = dict(name='RF', length=0.18856, voltage=30.0,
          phase=-1.5708+0.3491, frequency=704e6,
          cavity_radius=0.163, cavity_thickness=0.003, cavity_material='G4_Cu')

channel = build_coil_beamline(z_coils=[0.081, 0.211],
                              coil_templates=[coil_2, coil_1],
                              polarities=[1, 1], fixed=True, channel=channel)
channel = build_dipole_beamline(coil_cell_z=0.025, dipole_template=dipole, channel=channel)
channel = build_absorber_beamline(absorber_template=absorber, channel=channel)
channel = build_rf_beamline(n_rf_cells=3, rf_spacing=0.1946, rf_template=rf, channel=channel)

channel.summary()

render_gmad(channel, "config/channel.tpl", "output/channel.gmad",
            sampler_mode="linspace",
            sampler_kwargs={"n": 120, "start_m": -3.0, "end_m": 3.0},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim.dat"})
```

```bash
bdsim --file=output/channel.gmad --outfile=output --ngenerate=1000
```

See `quickstart.py` for the full 151-cell example.

---

## Physics

### Element classes

All elements share `z_center` (axial position [m]) and `name`. `fixed=True` locks the element during matching/optimisation.

#### SolenoidCoil

Models a solenoid as stacked cylindrical current sheets. Supports single-block or pancake-array geometry.

```python
SolenoidCoil(
    z_center    : float,
    r_in        : float,          # inner radius [m]
    r_thick     : float,          # radial thickness [m]

    # geometry — one of:
    L           : float = None,   # single block length [m]
    N_pancakes  : int   = 1,      # number of pancakes
    L_pancake   : float = None,   # pancake length [m]
    L_spacing   : float = 0.0,    # gap between pancakes [m]

    # current — one of:
    currDensity : float = None,   # [A/mm²]
    Itot        : float = None,   # total current [A]

    N_sheets    : int   = 10,     # radial integration points
    fixed       : bool  = False,
    material    : str   = "G4_Cu",
    name        : str   = None,
)
```

`get_field(z)` computes on-axis Bz using the exact current-sheet formula across all pancakes and radial sheets. Pass a negative current to reverse polarity.

#### Dipole

```python
Dipole(
    z_center         : float,
    field_strength   : float,
    aperture         : float = 0.2,
    length_z         : float = 0.1,
    enge_coefficient : float = 5.5,
    fixed            : bool  = False,
    name             : str   = None,
)
```

`get_field(z)` returns zero — the transverse dipole field is handled by BDSIM using the Enge fringe model.

#### Absorber

```python
Absorber(
    z_center              : float,
    absorber_type         : str   = "wedge",          # "wedge" or "cylinder"
    material              : str   = "G4_LITHIUM_HYDRIDE",
    cylinder_length       : float = 1.10,
    cylinder_radius       : float = 0.24,
    wedge_opening_angle   : float = 0.1745,           # [rad]
    wedge_height          : float = 0.24,
    wedge_rotation_angle  : float = 0.0,              # [rad] rotation about z
    wedge_offset_x        : float = 0.0,
    wedge_offset_y        : float = 0.0,
    wedge_apex_to_base    : float = 0.1134,
    fixed                 : bool  = False,
    name                  : str   = None,
)
```

`wedge_rotation_angle` alternates between cells (e.g. `[π, 0]`) to flip the wedge orientation and achieve alternating dispersion.

#### RFCavity

```python
RFCavity(
    z_center              : float,
    time_offset           : float = 0.0,              # transit time [ns] — set automatically
    length                : float = 0.18856,
    voltage               : float = 30.0,             # peak voltage [MV]
    phase                 : float = -π/2 + 0.349,
    frequency             : float = 704e6,
    window_thickness      : float = 0.0,
    window_material       : str   = "G4_Be",
    window_radius         : float = 0.0816,
    cavity_material       : str   = "G4_Cu",
    cavity_vacuum_material: str   = "vacuum",
    cavity_radius         : float = 0.163,
    cavity_thickness      : float = 0.003,
    fixed                 : bool  = False,
    name                  : str   = None,
)
```

`time_offset` is normally set automatically by `compute_rf_time_offsets()`.

---

### Beamline and MuonCoolingChannel

#### Beamline

Base container for an ordered sequence of elements.

```python
bl = Beamline(name="MyLine")
bl.add_element(elem)
bl.add_elements([e1, e2, e3])
bl.remove_element(index)
bl.clear()

bz     = bl.get_field(z_array)          # total Bz [T]
fields = bl.get_element_fields(z)       # dict of per-element fields
fig    = bl.plot_field(z_range=(a, b))
bl.summary()
```

#### MuonCoolingChannel

Extends `Beamline` with channel-level geometry and BDSIM field model settings.

```python
channel = MuonCoolingChannel(
    n_cells              : int,
    cell_length          : float,
    total_length         : float,
    total_width          : float,
    name                 : str   = "MuonCoolingChannel",
    magnetic_field_model : str   = "solenoidsheet",
    magnetic_field_method: str   = "grid",
    dipole_field_model   : str   = "dipole",
    interpolator         : str   = "linear",
    electric_field_model : str   = "rfpillbox",
    reference_momentum   : float = 200.0,             # [MeV/c]
    on_axis_tolerance    : float = 2e-2,
    polarities           : list[int] = None,          # default [1, -1, -1, 1]
)
```

**`compute_rf_time_offsets`** — synchronises each `RFCavity.time_offset` with a reference particle:

```python
offsets = channel.compute_rf_time_offsets(
    momentum_MeV_c = 200.0,
    beam_start     = 0.0,     # global z entry point [m]
    mass_MeV_c     = 105.66,
)
# returns {cavity_name: time_offset_ns, ...}
```

**`getCellStarts`** — returns `[[z_global, sign], ...]` at the start of each cell, with polarity sign derived from `self.polarities`.

---

### Builder functions

Each builder populates a `MuonCoolingChannel` with a periodically repeated element type. Pass `channel=` to append to an existing channel.

#### build_coil_beamline

```python
channel = build_coil_beamline(
    z_coils        : list[float],   # z positions from cell start [m]
    coil_templates : list[dict],    # one dict per coil type
    polarities     : list[int],     # +1 or -1 per coil type
    fixed          : bool = True,
    channel        : MuonCoolingChannel = None,
)
```

Each cell places two coils per entry: one at `z_local` and a mirror at `L_cell - z_local` with flipped polarity.

Template keys: `name`, `r_in`, `r_thick`, `N_pancakes`, `L_pancake`, `L_spacing`, `currDensity`, `N_sheets`, `material`

#### build_dipole_beamline

```python
channel = build_dipole_beamline(
    coil_cell_z     : float,
    dipole_template : dict,
    polarities      : list[int] = [1, -1, -1, 1],
    fixed           : bool = True,
    channel         : MuonCoolingChannel = None,
)
```

Template keys: `name`, `field_strength`, `aperture`, `length_z`, `enge_coefficient`

#### build_absorber_beamline

```python
channel = build_absorber_beamline(
    absorber_template : dict,
    rotation_angles   : list[float] = [π, 0.0],   # cycled per cell
    fixed             : bool = True,
    channel           : MuonCoolingChannel = None,
)
```

Places one absorber per cell at the cell centre. `wedge_rotation_angle` cycles through `rotation_angles`.

Template keys: `name`, `absorber_type`, `material`, `cylinder_length`, `cylinder_radius`, `wedge_opening_angle`, `wedge_height`, `wedge_apex_to_base`

#### build_rf_beamline

```python
channel = build_rf_beamline(
    n_rf_cells  : int,
    rf_spacing  : float,   # spacing between cavities [m]
    rf_template : dict,
    fixed       : bool = True,
    channel     : MuonCoolingChannel = None,
)
```

Places `n_rf_cells` cavities symmetrically about each cell centre.

Template keys: `name`, `length`, `voltage`, `phase`, `frequency`, `window_thickness`, `window_material`, `window_radius`, `cavity_material`, `cavity_vacuum_material`, `cavity_radius`, `cavity_thickness`

#### Low-level helpers

```python
positions = build_periodic_channel(N_cells, L_cell, z_coils, polarities)
z_list    = coil_placement_z(n_cells, cell_length, coil_cell_z)
```

---

### DataFrame introspection

All per-element DataFrames accept an optional `global_z_offset` keyword (default `0.0`).

```python
channel.build_pancake_dataframe()
# one row per sub-pancake: coil_name, coil_index, pancake_index,
# z_center, Itot_pancake, currDensity, polarity, r_in, r_thick, L_pancake, ...

channel.build_dipole_dataframe(global_z_offset=0.0)
# name, z_local, z_global, field_strength, aperture, length_z, enge_coefficient, fixed

channel.build_absorber_dataframe(global_z_offset=0.0)
# name, z_local, z_global, absorber_type, material, cylinder_*, wedge_*, fixed

channel.build_rf_dataframe(global_z_offset=0.0)
# name, z_local, z_global, time_offset, length, voltage, phase, frequency,
# window_*, cavity_*, fixed

channel.build_elements_dataframe(global_z_offset=0.0)
# name, type, z_local, z_global  (all element types, one row each)
```

---

## Rendering

### render_gmad

Renders a fully populated BDSIM `.gmad` from a `MuonCoolingChannel` and a Jinja2 template. Automatically syncs RF time offsets, builds the pancake DataFrame, formats all arrays as BDSIM `{...}` literals, and writes the file.

```python
render_gmad(
    beamline      : MuonCoolingChannel,
    tpl_path      : str,              # e.g. "config/channel.tpl"
    out_gmad      : str,              # e.g. "output/channel.gmad"
    sampler_mode  : str  = "linspace",  # "linspace" or "positions"
    sampler_kwargs: dict = None,
    beam_mode     : str  = "beam",    # "reference", "offset", or "beam"
    beam_kwargs   : dict = None,      # forwarded to generate_beam_block
)
```

**`sampler_mode="linspace"`** — evenly spaced samplers:

```python
render_gmad(channel, "config/channel.tpl", "output/channel.gmad",
            sampler_mode="linspace",
            sampler_kwargs={"n": 120, "start_m": -3.0, "end_m": 3.0,
                            "aper_m": 5.0},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim.dat"})
```

**`sampler_mode="positions"`** — arbitrary position array:

```python
import numpy as np
render_gmad(channel, "config/channel.tpl", "output/channel.gmad",
            sampler_mode="positions",
            sampler_kwargs={"positions_m": np.linspace(-3.0, 3.0, 120),
                            "aper_m": 5.0},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim.dat"})
```

`sampler_kwargs` keys for both modes:

| key | modes | default |
|---|---|---|
| `n` | `linspace` | 129 |
| `start_m` / `end_m` | `linspace` | ±`total_length/2` |
| `positions_m` | `positions` | required |
| `aper_m` | both | 5.0 |
| `reference_element` | both | `"mc1"` |
| `reference_element_number` | both | 0 |

The template (`config/channel.tpl`) is a Jinja2 file that defines the BDSIM material block, `cooldef1` (coolingchannel element), `mc1` (muoncooler), lattice line, options, and the injected `{{ sampler_lines }}` / `{{ beam_block }}` sections.

### generate_sampler_lines

Evenly spaced samplers via `numpy.linspace`:

```python
lines = generate_sampler_lines(
    n_samplers               : int,
    total_length_m           : float = None,
    aper_m                   : float = 5.0,
    reference_element        : str   = "mc1",
    reference_element_number : int   = 0,
    start_m                  : float = None,
    end_m                    : float = None,
)
# returns multi-line str of BDSIM samplerplacement commands
```

`start_m`/`end_m` take precedence over `total_length_m` when both are given.

### generate_sampler_lines_from_positions

Samplers at an arbitrary array of positions:

```python
lines = generate_sampler_lines_from_positions(
    positions_m              : array-like,  # positions [m]
    aper_m                   : float = 5.0,
    reference_element        : str   = "mc1",
    reference_element_number : int   = 0,
)
# returns multi-line str of BDSIM samplerplacement commands
```

### generate_beam_block

```python
block = generate_beam_block(
    mode             : str,          # "reference", "offset", or "beam"
    particle         : str   = "mu+",
    momentum_MeV     : float = 200.0,
    x0_mm            : float = 0.0,  # offset mode
    y0_mm            : float = 0.0,
    z0_m             : float = 0.0,
    distr_type       : str   = "userfile",
    distr_file       : str   = "",
    distr_file_format: str   = "t[ns]:x[m]:y[m]:z[m]:xp[rad]:yp[rad]:zp[rad]:-:E[GeV]",
    nlines_ignore    : int   = 1,
)
# returns str BDSIM beam command block
```

---

## Beam file translation

Converts G4Beamline `BLTrackFile` output to the `highemittanceBeamGen` / BDSIM `userfile` format.

```python
from src import g4bl_to_beamgen

n = g4bl_to_beamgen(
    input_path  : str,              # G4BL BLTrackFile (e.g. beam.tmp)
    output_path : str,
    mass_MeV    : float = 105.66,
    z_override  : float = None,     # pin all Z values to a fixed entry point [m]
)
# returns int — number of particles written
```

Input columns (G4BL): `x y z Px Py Pz t PDGid EvNum TrkId Parent weight` — positions in mm, momenta in MeV/c, time in ns.

Output columns: `T X Y Z xp yp zp P E` — T in ns, positions in m, direction cosines, P in GeV/c, E in GeV.

---

## JSON-driven pipeline

`src/pipeline.py` and `src/run.py` are a config-driven alternative to scripting the workflow by hand. Everything is described in `config/channel.json`; set any optional section to `null` to skip it.

```json
{
  "channel": {
    "n_cells": 151, "cell_length": 1.0, "total_length": 170.0,
    "total_width": 0.8, "reference_momentum": 200.0, "on_axis_tolerance": 2e-5
  },
  "coils": [
    {
      "name": "Coil2", "z_position": 0.081, "polarity": 1, "fixed": true,
      "r_in": 0.185, "r_thick": 0.060, "N_pancakes": 5,
      "L_pancake": 0.012, "L_spacing": 0.003, "currDensity": 300.0, "N_sheets": 5
    },
    {
      "name": "Coil1", "z_position": 0.211, "polarity": 1, "fixed": true,
      "r_in": 0.285, "r_thick": 0.070, "N_pancakes": 17,
      "L_pancake": 0.012, "L_spacing": 0.004, "currDensity": 328.43, "N_sheets": 5
    }
  ],
  "dipole": {
    "coil_cell_z": 0.025, "name": "Dipole", "field_strength": 0.2,
    "aperture": 0.2, "length_z": 0.1, "enge_coefficient": 5.5
  },
  "absorber": null,
  "rf": null,
  "output": {
    "template_path": "config/channel.tpl", "gmad_path": "output/channel.gmad",
    "output_name": "channel", "n_events": 100
  },
  "samplers": {
    "n_samplers": 120, "sampler_start_m": -3, "sampler_end_m": 3,
    "n_particles_per_sampler": 3
  },
  "beam": { "mode": "beam", "distr_file": "beam_bdsim.dat" },
  "analysis": { "output_csv": "output/sampler_data.csv" }
}
```

`coils` are grouped by `fixed` value — coils sharing the same `fixed` should be contiguous. `rf` must also include `n_rf_cells` and `rf_spacing`.

### pipeline.py — build only

```python
from src.pipeline import build_channel_from_config
import json

with open("config/channel.json") as f:
    config = json.load(f)

build_channel_from_config(config)
```

```bash
python -m src.pipeline --config config/channel.json
```

### run.py — build, simulate, and analyse

Calls `build_channel_from_config`, runs BDSIM via `pybdsim.Run.Bdsim`, then extracts sampler data to CSV.

```python
from src.run import run_from_config
import json

with open("config/channel.json") as f:
    config = json.load(f)

run_from_config(config)
```

```bash
python -m src.run --config config/channel.json
```
