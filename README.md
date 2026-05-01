# muon-cooling

Python toolkit for designing and rendering muon ionization cooling channels as BDSIM input files. It provides a class hierarchy for beamline elements, periodic-channel builder functions, and a Jinja2-based GMAD renderer.

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Element Classes](#element-classes)
   - [SolenoidCoil](#solenoidcoil)
   - [Dipole](#dipole)
   - [Absorber](#absorber)
   - [RFCavity](#rfcavity)
4. [Beamline and MuonCoolingChannel](#beamline-and-muoncoolingchannel)
5. [Builder Functions](#builder-functions)
6. [DataFrame Introspection](#dataframe-introspection)
7. [GMAD Rendering](#gmad-rendering)
8. [Template System](#template-system)
9. [Quickstart](#quickstart)

---

## Overview

A muon ionization cooling channel consists of repeating cells, each containing:

- **Solenoid coils** — provide the focusing magnetic field
- **Wedge/cylinder absorbers** — reduce muon momentum (cooling)
- **RF cavities** — restore longitudinal momentum
- **Dipoles** — provide dispersion for emittance exchange

This library models these elements in Python, lets you build periodic channels programmatically, and renders them to BDSIM-compatible `.gmad` files via a Jinja2 template.

---

## Architecture

```
BeamlineElement          (base class)
├── SolenoidCoil         (axial magnetic field from stacked current sheets)
├── Dipole               (transverse dipole field)
├── Absorber             (material absorber: wedge or cylinder)
└── RFCavity             (RF accelerating cavity)

Beamline                 (ordered list of BeamlineElement objects)
└── MuonCoolingChannel   (Beamline + channel-level geometry and field model settings)
```

All elements share a `z_center` (axial position in meters) and a `name`. The `Beamline` class superimposes their fields along the z axis.

---

## Element Classes

### SolenoidCoil

Models a solenoid as stacked cylindrical current sheets. Supports single-block or pancake-array geometry.

```python
SolenoidCoil(
    z_center      : float,          # axial center [m]
    r_in          : float,          # inner radius [m]
    r_thick       : float,          # radial thickness [m]

    # Geometry — choose one:
    L             : float = None,   # single block length [m]
    N_pancakes    : int   = 1,      # number of pancakes
    L_pancake     : float = None,   # pancake length [m]
    L_spacing     : float = 0.0,    # gap between pancakes [m]

    # Current — choose one:
    currDensity   : float = None,   # current density [A/mm²]
    Itot          : float = None,   # total current [A]

    N_sheets      : int   = 10,     # radial integration points
    fixed         : bool  = False,  # hold current fixed during matching
    material      : str   = "G4_Cu",
    name          : str   = None,
)
```

**Field calculation** — `get_field(z)` sums the on-axis Bz contribution of each pancake and each radial sheet using the exact current-sheet formula:

```
Bz(z) = (μ₀ I) / (2 L) * [ (z+ / √(z+² + r²)) - (z- / √(z-² + r²)) ]
```

where `z± = z - z0 ± L/2`.

**Properties:**
- `Itot` — total current across all pancakes (read/write)
- `currDensity_value` — current density [A/mm²] (read/write)

**Polarity** — pass a negative `currDensity` or `Itot` to reverse the field direction.

---

### Dipole

A transverse dipole element. Field is represented only by its strength; the actual field map is handled by BDSIM using the Enge fringe model.

```python
Dipole(
    z_center         : float,       # axial center [m]
    field_strength   : float,       # dipole field [T]
    aperture         : float = 0.2, # aperture half-gap [m]
    length_z         : float = 0.1, # longitudinal length [m]
    enge_coefficient : float = 5.5, # Enge fringe coefficient
    fixed            : bool  = False,
    name             : str   = None,
)
```

`get_field(z)` returns zero (Bz only; the dipole's transverse field is handled downstream by BDSIM).

---

### Absorber

A material absorber, either a **wedge** (for ionization cooling with emittance exchange) or a **cylinder** (for pure longitudinal cooling).

```python
Absorber(
    z_center              : float,
    absorber_type         : str   = "wedge",             # "wedge" or "cylinder"
    material              : str   = "G4_LITHIUM_HYDRIDE",

    # Cylinder geometry:
    cylinder_length       : float = 1.10,   # [m]
    cylinder_radius       : float = 0.24,   # [m]

    # Wedge geometry:
    wedge_opening_angle   : float = 0.1745, # [rad] (≈10°)
    wedge_height          : float = 0.24,   # [m]
    wedge_rotation_angle  : float = 0.0,    # [rad] rotation about z
    wedge_offset_x        : float = 0.0,    # [m]
    wedge_offset_y        : float = 0.0,    # [m]
    wedge_apex_to_base    : float = 0.1134, # [m]

    fixed                 : bool  = False,
    name                  : str   = None,
)
```

`wedge_rotation_angle` alternates between cells (e.g. `[π, 0]`) to flip the wedge orientation and achieve alternating dispersion.

---

### RFCavity

An RF pillbox cavity that restores longitudinal momentum lost in the absorbers.

```python
RFCavity(
    z_center              : float,
    time_offset           : float = 0.0,      # transit time offset [ns]
    length                : float = 0.18856,  # cavity length [m]
    voltage               : float = 30.0,     # peak voltage [MV]
    phase                 : float = -π/2 + 0.349, # RF phase [rad]
    frequency             : float = 704e6,    # RF frequency [Hz]
    window_thickness      : float = 0.0,      # entrance window thickness [m]
    window_material       : str   = "G4_Be",
    window_radius         : float = 0.0816,   # [m]
    cavity_material       : str   = "G4_Cu",
    cavity_vacuum_material: str   = "vacuum",
    cavity_radius         : float = 0.163,    # [m]
    cavity_thickness      : float = 0.003,    # [m]
    fixed                 : bool  = False,
    name                  : str   = None,
)
```

`time_offset` is normally set automatically by `Beamline.compute_rf_time_offsets()` — see below.

---

## Beamline and MuonCoolingChannel

### Beamline

The base container for an ordered sequence of elements.

```python
bl = Beamline(name="MyLine")
bl.add_element(elem)
bl.add_elements([e1, e2, e3])
bl.remove_element(index)
bl.clear()

bz = bl.get_field(z_array)            # total Bz [T] at each z
fields = bl.get_element_fields(z)     # dict of per-element fields
fig, ax = bl.plot_field(z_range=(z_min, z_max), n_points=1000)
bl.summary()                          # print element table
```

#### compute_rf_time_offsets

Sets `time_offset` on every `RFCavity` to synchronize it with a reference particle:

```python
offsets = bl.compute_rf_time_offsets(
    momentum_MeV_c = 200.0,   # reference muon momentum [MeV/c]
    mass_MeV_c     = 105.66,  # muon rest mass [MeV/c²]
    z_ref          = None,    # reference z; defaults to leftmost element
)
# returns dict: {cavity_name: time_offset_ns, ...}
```

Transit time is computed as `t = (z_cavity - z_ref) / (β c)`.

---

### MuonCoolingChannel

Extends `Beamline` with channel-level geometry and BDSIM field model settings.

```python
channel = MuonCoolingChannel(
    n_cells              : int,
    cell_length          : float,   # [m]
    total_length         : float,   # full channel length [m]
    total_width          : float,   # horizontal width [m]
    name                 : str   = "MuonCoolingChannel",
    magnetic_field_model : str   = "solenoidsheet",
    magnetic_field_method: str   = "grid",
    dipole_field_model   : str   = "dipole",
    interpolator         : str   = "linear",
    electric_field_model : str   = "rfpillbox",
    reference_momentum   : float = 200.0,  # [MeV/c]
    on_axis_tolerance    : float = 2e-2,
)
```

---

## Builder Functions

These functions populate a `MuonCoolingChannel` with elements arranged in a periodic pattern.

### build_coil_beamline — solenoid coils

```python
channel = build_coil_beamline(
    z_coils        : list[float],   # coil z positions within a cell from cell start [m]
    coil_templates : list[dict],    # one dict per coil type (see below)
    polarities     : list[int],     # +1 or -1 per coil type
    fixed          : bool = True,
    channel        : MuonCoolingChannel = None,  # append to existing channel
)
```

Each cell generates two coils per entry in `z_coils`: one at `z_local` and a mirror at `L_cell - z_local` with flipped polarity.

**Coil template keys:** `name`, `r_in`, `r_thick`, `N_pancakes`, `L_pancake`, `L_spacing`, `currDensity`, `N_sheets`, `material`

---

### build_dipole_beamline — dipoles

```python
channel = build_dipole_beamline(
    coil_cell_z     : float,        # dipole offset from cell edge [m]
    dipole_template : dict,         # see Dipole parameters
    polarities      : list[int] = [1, -1, -1, 1],
    fixed           : bool = True,
    channel         : MuonCoolingChannel = None,
)
```

**Dipole template keys:** `name`, `field_strength`, `aperture`, `length_z`, `enge_coefficient`

---

### build_absorber_beamline — absorbers

```python
channel = build_absorber_beamline(
    absorber_template : dict,
    rotation_angles   : list[float] = [π, 0.0],  # cycled per cell
    fixed             : bool = True,
    channel           : MuonCoolingChannel = None,
)
```

Places one absorber per cell at the cell center (`z_cell_start`). The `wedge_rotation_angle` cycles through `rotation_angles`.

**Absorber template keys:** `name`, `absorber_type`, `material`, `cylinder_length`, `cylinder_radius`, `wedge_opening_angle`, `wedge_height`, `wedge_apex_to_base`

---

### build_rf_beamline — RF cavities

```python
channel = build_rf_beamline(
    n_rf_cells  : int,      # cavities per cell
    rf_spacing  : float,    # spacing between cavities [m]
    rf_template : dict,
    fixed       : bool = True,
    channel     : MuonCoolingChannel = None,
)
```

Places `n_rf_cells` cavities symmetrically about each cell center, separated by `rf_spacing`.

**RF template keys:** `name`, `length`, `voltage`, `phase`, `frequency`, `window_thickness`, `window_material`, `window_radius`, `cavity_material`, `cavity_vacuum_material`, `cavity_radius`, `cavity_thickness`

---

### Low-level helpers

```python
# Raw list of {z, polarity, coil_type} dicts for a periodic layout:
positions = build_periodic_channel(N_cells, L_cell, z_coils, polarities)

# z positions for coil_placement_z layout (legacy helper):
z_list = coil_placement_z(n_cells, cell_length, coil_cell_z)
```

---

## DataFrame Introspection

Each element type has a corresponding DataFrame builder on `Beamline`/`MuonCoolingChannel`. These are useful for inspection, export, and debugging — the render pipeline uses `build_pancake_dataframe` internally.

### build_pancake_dataframe

One row per sub-pancake across all `SolenoidCoil` elements.

```python
df = channel.build_pancake_dataframe()
# Columns: coil_name, coil_index, pancake_index, z_center, z_coil_center,
#          Itot_pancake, currDensity, polarity, r_in, r_thick, L_pancake,
#          fixed, material, coil_offset_x, coil_offset_y,
#          coil_tilt_x, coil_tilt_y, coil_tilt_z
```

### build_dipole_dataframe

One row per `Dipole`.

```python
df = channel.build_dipole_dataframe()
# Columns: name, element_index, z_center, field_strength, aperture,
#          length_z, enge_coefficient, fixed
```

### build_absorber_dataframe

One row per `Absorber`.

```python
df = channel.build_absorber_dataframe()
# Columns: name, element_index, z_center, absorber_type, material,
#          cylinder_length, cylinder_radius, wedge_opening_angle,
#          wedge_height, wedge_rotation_angle, wedge_offset_x,
#          wedge_offset_y, wedge_apex_to_base, fixed
```

### build_rf_dataframe

One row per `RFCavity`.

```python
df = channel.build_rf_dataframe()
# Columns: name, element_index, z_center, time_offset, length, voltage,
#          phase, frequency, window_thickness, window_material,
#          window_radius, cavity_material, cavity_vacuum_material,
#          cavity_radius, cavity_thickness, fixed
```

---

## GMAD Rendering

`render_gmad` writes a fully populated BDSIM `.gmad` file from a `MuonCoolingChannel` and a Jinja2 template.

```python
render_gmad(
    beamline     : MuonCoolingChannel,
    tpl_path     : str,             # path to templates/channel.tpl
    out_gmad     : str,             # output path
    n_samplers   : int   = 129,     # number of samplerplacement elements
    sampler_aper_m: float = 5.0,    # sampler aperture half-size [m]
    beam_mode    : str   = "beam",  # "reference", "offset", or "beam"
    beam_kwargs  : dict  = None,    # extra kwargs forwarded to generate_beam_block
)
```

It automatically:
1. Calls `compute_rf_time_offsets` to synchronize RF phases.
2. Builds the pancake DataFrame and sorts it by coil/pancake index.
3. Collects dipoles, absorbers, and RF cavities from the element list.
4. Formats all arrays as BDSIM `{...}` literals (collapsed to a scalar when uniform).
5. Renders the Jinja2 template and writes the output file.

### generate_sampler_lines

```python
lines = generate_sampler_lines(
    n_samplers               : int,
    total_length_m           : float,
    aper_m                   : float = 5.0,
    reference_element        : str   = "mc1",
    reference_element_number : int   = 0,
)
```

Returns a multi-line string of `samplerplacement` BDSIM commands evenly spaced across `±total_length_m/2`.

### generate_beam_block

```python
block = generate_beam_block(
    mode            : str,          # "reference", "offset", or "beam"
    particle        : str   = "mu+",
    momentum_MeV    : float = 200.0,
    # offset mode:
    x0_mm, y0_mm   : float = 0.0,
    z0_m            : float = 0.0,
    # beam mode:
    distr_type      : str   = "userfile",
    distr_file      : str   = "",
    distr_file_format: str  = "t[ns]:x[m]:y[m]:z[m]:xp[rad]:yp[rad]:zp[rad]:-:E[GeV]",
    nlines_ignore   : int   = 1,
)
```

---

## Template System

`templates/channel.tpl` is a Jinja2 template with `{{ variable }}` placeholders. The template defines:

- **Material definitions** (custom BDSIM materials)
- **`cooldef1`** — the `coolingchannel` BDSIM element, parameterized by all coil/dipole/absorber/RF arrays
- **`mc1`** — the `muoncooler` element referencing `cooldef1`, with total length and width
- **Lattice line** — `lat: line=(mc1, d1, c1)`
- **Sampler placements** — injected via `{{ sampler_lines }}`
- **Options block** — physics list, step limits, overlap checking
- **Beam block** — injected via `{{ beam_block }}`

To customize the template, edit `templates/channel.tpl` directly. All scalar/array values passed to `render_gmad` are available in the template context.

---

## Quickstart

```python
from src import (
    MuonCoolingChannel,
    build_coil_beamline,
    build_dipole_beamline,
    build_absorber_beamline,
    build_rf_beamline,
    render_gmad,
)

# 1. Create the channel
channel = MuonCoolingChannel(
    n_cells            = 10,
    cell_length        = 0.8,
    total_length       = 124.0,
    total_width        = 0.8,
    reference_momentum = 200.0,
    on_axis_tolerance  = 2e-2,
)

# 2. Define element templates
coil_1_template = {
    'name': 'Coil1', 'r_in': 0.285, 'r_thick': 0.070,
    'N_pancakes': 17, 'L_pancake': 0.012, 'L_spacing': 0.004,
    'currDensity': 328.43, 'N_sheets': 5,
}
coil_2_template = {
    'name': 'Coil2', 'r_in': 0.185, 'r_thick': 0.060,
    'N_pancakes': 5, 'L_pancake': 0.012, 'L_spacing': 0.003,
    'currDensity': 300.0, 'N_sheets': 5,
}
dipole_template = {
    'name': 'Dipole', 'field_strength': 0.2,
    'aperture': 0.2, 'length_z': 0.1, 'enge_coefficient': 5.5,
}
absorber_template = {
    'name': 'Absorber', 'absorber_type': 'wedge',
    'material': 'G4_LITHIUM_HYDRIDE',
    'wedge_opening_angle': 0.1745, 'wedge_height': 0.24,
    'wedge_apex_to_base': 0.1134256,
}
rf_template = {
    'name': 'RF', 'length': 0.18856, 'voltage': 30.0,
    'phase': -1.5708 + 0.3491, 'frequency': 704e6,
    'cavity_radius': 0.163, 'cavity_thickness': 0.003,
    'cavity_material': 'G4_Cu',
}

# 3. Populate channel with elements
channel = build_coil_beamline(
    z_coils=[0.081, 0.211], coil_templates=[coil_2_template, coil_1_template],
    polarities=[1, 1], fixed=True, channel=channel,
)
channel = build_dipole_beamline(coil_cell_z=0.7, dipole_template=dipole_template, channel=channel)
channel = build_absorber_beamline(absorber_template=absorber_template, channel=channel)
channel = build_rf_beamline(n_rf_cells=3, rf_spacing=0.1946, rf_template=rf_template, channel=channel)

# 4. Inspect
channel.summary()
print(channel.build_pancake_dataframe())
print(channel.build_rf_dataframe())

# 5. Render to GMAD
render_gmad(channel, "templates/templates/channel.tpl", "channel.gmad", n_samplers=10)
```

This produces `channel.gmad`, which can be run directly with BDSIM:

```bash
bdsim --file=channel.gmad --outfile=output --ngenerate=1000
```
