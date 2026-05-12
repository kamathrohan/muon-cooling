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
8. [Condor job setup](#condor-job-setup)
9. [Analysis](#analysis)
   - [loader.py](#loaderpy)
   - [optics.py](#opticspy)
   - [bdsim_to_g4bl](#bdsim_to_g4bl)
   - [plot_channel.py](#plot_channelpy)

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
│   └── pipeline.py       # config-driven build (JSON → gmad)
├── analyse/
│   ├── loader.py         # file discovery, ROOT reading, CSV loading
│   ├── optics.py         # emittance, beta, optics calculations and plots
│   └── bdsim_to_g4bl.py  # convert BDSIM DataFrame to G4BL/ICOOL track format
├── condor/
│   ├── helper.py         # Condor script rendering and folder utilities
│   └── setup_sim.py      # job initialisation (folders, gmad, beam file)
├── config/               # Jinja2 templates and JSON configs
├── output/               # simulation artefacts (gitignored)
├── quickstart.py         # channel build example
└── plot_channel.py       # analysis CLI and importable run() function
```

---

## Quickstart

```python
import json
import numpy as np
from src import (
    MuonCoolingChannel,
    build_coil_beamline, build_dipole_beamline,
    build_absorber_beamline, build_rf_beamline,
    render_gmad,
)

channel = MuonCoolingChannel(
    n_cells=151, cell_length=1.0, total_length=170.0,
    total_width=0.8, reference_momentum=200.0, on_axis_tolerance=2e-5,
    magnetic_field_model="solenoidblock", magnetic_field_method="grid",
    grid_points_per_mm=1.0,
    z_period_start=-70000, z_period_end=175000, period_length=2000,
)

coil_2 = dict(name='Coil2', r_in=0.185, r_thick=0.060,
              N_pancakes=5, L_pancake=0.012, L_spacing=0.003,
              currDensity=300.0, N_sheets=5)
coil_1 = dict(name='Coil1', r_in=0.285, r_thick=0.070,
              N_pancakes=17, L_pancake=0.012, L_spacing=0.004,
              currDensity=328.43, N_sheets=5)
dipole   = dict(name='Dipole', field_strength=0.2,
                aperture=0.2, length_z=0.1, enge_coefficient=5.5)
absorber = dict(name='Absorber', absorber_type='wedge',
                material='G4_LITHIUM_HYDRIDE',
                wedge_opening_angle=0.1745, wedge_height=0.28571,
                wedge_apex_to_base=0.28571)
rf = dict(name='RF', length=0.18856, voltage=30.0,
          phase=-1.5707963267948966 + 0.3490658503988659, frequency=704e6,
          cavity_radius=0.163, cavity_thickness=0.003, cavity_material='G4_Cu')

channel = build_coil_beamline(z_coils=[0.081, 0.211],
                              coil_templates=[coil_2, coil_1],
                              polarities=[1, 1], fixed=True, channel=channel)
channel = build_dipole_beamline(coil_cell_z=0.025, dipole_template=dipole, channel=channel)
channel = build_absorber_beamline(n_cells=129, absorber_template=absorber,
                                  wedge_alignment_angle1=0.0, wedge_alignment_angle2=np.pi,
                                  channel=channel)
channel = build_rf_beamline(n_cells=129, n_rf_cells=3, rf_spacing=0.1946,
                            rf_template=rf, channel=channel)

with open("config/rfPhases.json") as f:
    rf_phases = json.load(f)
channel.set_rf_time_offsets(rf_phases["closedOrbit"])

render_gmad(channel, "config/channel.tpl", "output/channel.gmad",
            sampler_mode="linspace",
            sampler_kwargs={"n": 200, "start_m": -65, "end_m": 65},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim.dat"})
```

```bash
bdsim --file=output/channel.gmad --outfile=output --ngenerate=1000
```

See `quickstart.py` for the full working example.

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

    # misalignment
    tilt_x      : float = 0.0,   # [rad]
    tilt_y      : float = 0.0,   # [rad]
    tilt_z      : float = 0.0,   # [rad]
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

`wedge_rotation_angle` alternates between cells to flip the wedge orientation and achieve alternating dispersion.

#### RFCavity

```python
RFCavity(
    z_center              : float,
    time_offset           : float = None,             # transit time [ns] — set before rendering
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

`time_offset` must be set before rendering, either via `compute_rf_time_offsets()` or `set_rf_time_offsets()`.

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
    grid_points_per_mm   : float = 0.1,
    z_period_start       : float = -70000,            # [mm]
    z_period_end         : float = 175000,            # [mm]
    period_length        : float = 2000,              # [mm]
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

**`set_rf_time_offsets`** — manually set time offsets from an array (e.g. loaded from a JSON file):

```python
channel.set_rf_time_offsets(offsets)   # array-like, one value per RFCavity in order
```

**`set_tilts`** — apply per-coil tilt angles to all `SolenoidCoil` elements:

```python
channel.set_tilts({
    "tiltX": [...],   # array of length n_coils [rad]
    "tiltY": [...],
    "tiltZ": [...],
})
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
    absorber_template      : dict,
    n_cells                : int   = None,        # defaults to channel.n_cells
    wedge_alignment_angle1 : float = math.pi,     # rotation for even sub-cells [rad]
    wedge_alignment_angle2 : float = 0.0,         # rotation for odd sub-cells [rad]
    offsetX                : float = 0.0,         # x-offset magnitude; sign alternates [m]
    fixed                  : bool  = True,
    channel                : MuonCoolingChannel = None,
)
```

Places one absorber per cell at the cell centre. Wedge rotation strictly alternates by cell index parity.

Template keys: `name`, `absorber_type`, `material`, `cylinder_length`, `cylinder_radius`, `wedge_opening_angle`, `wedge_height`, `wedge_apex_to_base`

#### build_rf_beamline

```python
channel = build_rf_beamline(
    n_rf_cells  : int,
    rf_spacing  : float,            # spacing between cavities [m]
    rf_template : dict,
    n_cells     : int   = None,     # defaults to channel.n_cells
    fixed       : bool  = True,
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
# z_center, z_coil_center, Itot_pancake, currDensity, polarity,
# r_in, r_thick, L_pancake, fixed, material,
# coil_offset_x, coil_offset_y, coil_tilt_x, coil_tilt_y, coil_tilt_z

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

Renders a fully populated BDSIM `.gmad` from a `MuonCoolingChannel` and a Jinja2 template. Validates element bounds and that all RF cavities have `time_offset` set, builds the pancake DataFrame, formats all arrays as BDSIM `{...}` literals, and writes the file.

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
            sampler_kwargs={"n": 200, "start_m": -65, "end_m": 65,
                            "aper_m": 5.0},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim.dat"})
```

**`sampler_mode="positions"`** — arbitrary position array:

```python
import numpy as np
render_gmad(channel, "config/channel.tpl", "output/channel.gmad",
            sampler_mode="positions",
            sampler_kwargs={"positions_m": np.linspace(-65, 65, 200),
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

`src/pipeline.py` is a config-driven alternative to scripting the workflow by hand. Everything is described in `config/channel.json`; set any optional section to `null` to skip it.

```json
{
  "channel": {
    "n_cells": 151, "cell_length": 1.0, "total_length": 170.0,
    "total_width": 0.8, "reference_momentum": 200.0, "on_axis_tolerance": 2e-5,
    "magnetic_field_model": "solenoidblock", "magnetic_field_method": "grid",
    "grid_points_per_mm": 1.0,
    "z_period_start": -70000, "z_period_end": 175000, "period_length": 2000
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
  "absorber": {
    "name": "Absorber", "absorber_type": "wedge", "material": "G4_LITHIUM_HYDRIDE",
    "wedge_opening_angle": 0.1745, "wedge_height": 0.28571, "wedge_apex_to_base": 0.28571,
    "wedge_alignment_angle1": 0.0, "wedge_alignment_angle2": 3.14159,
    "wedge_offset_x": 0.0, "n_cells": 129
  },
  "rf": {
    "name": "RF", "length": 0.18856, "voltage": 30.0, "frequency": 704e6,
    "cavity_radius": 0.163, "cavity_thickness": 0.003, "cavity_material": "G4_Cu",
    "n_rf_cells": 3, "rf_spacing": 0.1946, "n_cells": 129
  },
  "rf_phases": [0.0, 1.23, ...],
  "coilTilts": {
    "tiltX": 0.0, "tiltY": 0.0, "tiltZ": 0.0
  },
  "tolerance": {
    "coil": { "current": 0.001, "tilt": 0.0005, "offset": 0.0, "seed": 42 }
  },
  "samplers": {
    "sampler_mode": "linspace",
    "sampler_kwargs": { "n": 200, "start_m": -65, "end_m": 65 }
  },
  "beam": { "mode": "beam", "distr_file": "beam_bdsim.dat" }
}
```

`coils` are grouped by `fixed` value — coils sharing the same `fixed` should be contiguous. `rf` must include `n_rf_cells` and `rf_spacing`. `rf_phases` and `rf_phasing` are mutually exclusive: use `rf_phases` for an explicit array of time offsets, or `rf_phasing: {"mode": "compute", "momentum_MeV_c": ..., "beam_start": ...}` to compute them automatically. `coilTilts` accepts either a scalar (applied to all coils) or an array of length equal to the number of solenoid coils.

### pipeline.py — build only

```python
from src.pipeline import build_channel_from_config
import json

with open("config/channel.json") as f:
    config = json.load(f)

build_channel_from_config(config, "config/channel.tpl", "output/channel.gmad")
```

```bash
python -m src.pipeline --config config/channel.json \
    --template config/channel.tpl --output output/channel.gmad
```

---

## Condor job setup

`condor/helper.py` and `condor/setup_sim.py` automate setting up batches of HTCondor simulation jobs.

### setup_simulation

Renders a single simulation directory with a `.gmad` file and beam distribution copy:

```python
from condor.setup_sim import setup_simulation

gmad_path = setup_simulation(
    simDir              : str,   # directory to create
    channel_json        : str,   # path to channel JSON config
    template            : str,   # path to Jinja2 .tpl file
    beam_file           : str,   # beam distribution to copy in
    gmad_name           : str  = None,    # override output filename
    distr_file_override : str  = None,    # override path in gmad (skip copy)
    extra_config        : dict = None,    # shallow-merged into config before render
)
# returns path to rendered .gmad
```

### initialiseJob

Sets up all batch/replica folders, renders the `.gmad` into each one, and writes the Condor submission scripts:

```python
from condor.setup_sim import initialiseJob

replica_folders = initialiseJob(
    simDir, outputDir, jobId, iterNum,
    nBatch, nReplicas,
    channel_json, beam_file, channel_template,
    datagen_template, submit_template, do_datagen_template,
    max_runtime, env_setup, bdsim_setup, ngenerate,
    tolerance = {"coil": {"current": 0.001, "tilt": 0.0, "offset": 0.0}},
    split     = False,   # if True, divides ngenerate across replicas
)
```

Each batch draws a fresh random seed for coil tolerances. `split=True` divides `ngenerate` evenly across replicas and staggers `skiplines` so each replica reads a disjoint slice of the beam file.

```bash
python condor/setup_sim.py --config config/job.json --job-id run01 [--split] [--run]
```

`--run` immediately executes `doDataGeneration.sh` after setup.

---

## Analysis

Post-processing tools live in `analyse/`. They are independent of the `src/` build pipeline and work on BDSIM sampler output.

### loader.py

File discovery, ROOT reading, and CSV loading.

```python
from analyse.loader import find_files, load_sims, analyze_root_file

# find files
txt_files = find_files("/path/to/sims")               # default pattern "*.txt"
csv_files = find_files("/path/to/sims", "*.csv")

# load and concatenate — tags each row with SimName from the filename stem
df = load_sims(txt_files)

# read a BDSIM ROOT file into a DataFrame (saves .csv alongside by default)
df, csv_path = analyze_root_file("output/channel.root")
df           = analyze_root_file("output/channel.root", save=False)
```

`load_sims` expects each file to be a CSV with at least `S`, `x`, `y`, `z`, `xp`, `yp`, `zp`, `p`, `energy`, `T`, `weight`, `trackID`, `samplerName` columns.

### optics.py

Emittance, beta function, and optics calculations, plus publication-quality plots.

#### Cuts

```python
from analyse.optics import apply_cuts

df = apply_cuts(
    df,
    r_max       = 0.0819,    # transverse aperture cut [m]; None to disable
    s_max       = 100000,    # keep S <= s_max [mm]; None to disable
    transmitted = True,      # keep only particles reaching the last station
)
```

#### Optics calculation

```python
from analyse.optics import compute_optics

optics = compute_optics(df)
# columns: samplerName, mean_z [mm], emittance_trans, emittance_long, emittance_6d,
#          beta_x, beta_y [mm], mean_p [MeV/c], transmission [%], a_pz_corr
```

All emittances are normalised RMS values. Units: transverse/longitudinal in mm·rad, 6D in mm³·rad³.

#### Individual physics functions

```python
from analyse.optics import (
    transverse_emittance,    # normalised 4D transverse emittance [mm]
    longitudinal_emittance,  # normalised longitudinal emittance [mm]
    emittance_6d,            # normalised 6D emittance [mm³]
    beta_function,           # (beta_x, beta_y) [mm]
    amplitude_pz_correlation,
)

eps_t = transverse_emittance(group)     # group is a per-station DataFrame slice
beta_x, beta_y = beta_function(group)
```

#### Plotting

```python
from analyse.optics import (
    plot_beta, plot_emittance, plot_transmission,
    plot_momentum, plot_amplitude_pz_corr,
    plot_all, make_phase_space_movie,
)

plot_beta(optics, save_path="beta.png")
plot_emittance(optics, save_path="emittance.png")
plot_transmission(df_raw, r_max=0.0819, save_path="transmission.png")

# save all five figures at once
plot_all(optics, df_raw, show=False, save="plots/")

# phase-space GIF (one frame per station)
make_phase_space_movie(df, out_path="plots/phase_space.gif", fps=5)
```

### bdsim_to_g4bl

Converts a BDSIM-style DataFrame to G4BL/ICOOL track format for use with xboa.

```python
from analyse.bdsim_to_g4bl import bdsim_to_g4bl

station_map = bdsim_to_g4bl(df, "output/output.txt")
```

Input DataFrame must have: `samplerName`, `trackID`, `SimName`, `partID`, `x`, `y`, `z`, `xp`, `yp`, `zp`, `p`, `T`, `weight`.

Output format: space-separated, one particle per line, with ICOOL PID codes and time converted to seconds.

### plot_channel.py

Top-level script that wires loader → optics → plots. Can be used from the CLI or imported into a notebook.

```python
from plot_channel import run

optics, df, dfAll = run(
    "/path/to/sims",
    pattern     = "*.txt",
    s_offset    = 21.5,      # subtracted from S [mm]
    r_max       = 0.0819,    # aperture cut [m]
    s_max       = 100000,    # zoom to channel region [mm]
    transmitted = True,
    show        = True,      # display interactively
    save        = "plots/",  # save to disk; None to skip
    movie       = False,
)
```

```bash
python plot_channel.py /path/to/sims --s-offset 21.5 --r-max 0.0819 --save plots/
python plot_channel.py /path/to/sims --show --no-save
python plot_channel.py /path/to/sims --movie
```
