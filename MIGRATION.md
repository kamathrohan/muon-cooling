# CoolingPancakes → src/ migration cheatsheet

Mechanical find-and-replace guide for any script importing from `CoolingPancakes`.

---

## 1. Import line

```python
# BEFORE
from CoolingPancakes import SolenoidCoil, Beamline, build_periodic_channel, build_coil_beamline

# AFTER
from src import (
    SolenoidCoil, Beamline, MuonCoolingChannel,
    build_periodic_channel, build_coil_beamline,
    build_dipole_beamline, build_absorber_beamline, build_rf_beamline,
    render_gmad,
)
```

---

## 2. Replace `Beamline` with `MuonCoolingChannel`

Wherever you created a top-level `Beamline` to hold the full cooling channel,
replace it with `MuonCoolingChannel`. `MuonCoolingChannel` is a subclass of
`Beamline` so all existing method calls (`add_element`, `get_field`,
`plot_field`, `summary`, etc.) still work unchanged.

```python
# BEFORE
beamline = Beamline(name="7-Cell Cooling Channel")

# AFTER
channel = MuonCoolingChannel(
    n_cells            = 7,
    cell_length        = 1.0,
    total_length       = 7.0,      # full channel length [m]
    total_width        = 0.8,      # [m]
    reference_momentum = 200.0,    # [MeV/c]
    on_axis_tolerance  = 2e-2,
)
```

---

## 3. `build_channel_beamline` renamed to `build_coil_beamline` — call pattern also changed

The old version took positional arguments and returned a plain `Beamline`.
The new version is renamed, takes keyword arguments, accepts an existing `MuonCoolingChannel`
via `channel=`, and returns that same channel.

```python
# BEFORE
beamline = build_channel_beamline(
    N_cells,
    L_cell,
    z_coils,
    coil_templates,
    polarities=[1, 1],
    fixed=True,
    name="My Channel",
)

# AFTER  (create the channel first, then populate it)
channel = MuonCoolingChannel(n_cells=N_cells, cell_length=L_cell, ...)
channel = build_coil_beamline(
    z_coils        = z_coils,
    coil_templates = coil_templates,
    polarities     = [1, 1],
    fixed          = True,
    channel        = channel,   # N_cells and L_cell are read from here
)
```

---

## 4. Coil templates — add `material` if overriding

The `SolenoidCoil` constructor gained a `material` parameter (default
`"G4_Cu"`). Existing template dicts and `SolenoidCoil(...)` calls work
unchanged; only add `"material"` if you need a different conductor material.

---

## 5. Summary of name changes

| Old | New |
|---|---|
| `from CoolingPancakes import ...` | `from src import ...` |
| `Beamline(name=...)` as top-level container | `MuonCoolingChannel(n_cells, cell_length, total_length, total_width, reference_momentum, ...)` |
| `build_channel_beamline(N, L, z, tmpls, ...)` positional | `build_coil_beamline(z_coils=..., coil_templates=..., channel=channel)` kwargs — **renamed + new signature** |
| No GMAD output | `render_gmad(channel, "templates/channel.tpl", "channel.gmad")` |
