import math
import numpy as np

from .elements import SolenoidCoil, Dipole, Absorber, RFCavity
from .beamline import MuonCoolingChannel


def coil_placement_z(n_cells, cell_length, coil_cell_z):
    coil_global_z = []
    if n_cells % 2 != 0:
        for cell_index in range(n_cells):
            z_0 = (-n_cells/2 + 0.5 + cell_index) * cell_length
            z_coil1 = z_0 - cell_length/2 + coil_cell_z
            z_coil2 = z_0 + cell_length/2 - coil_cell_z
            coil_global_z.append(z_coil1)
            coil_global_z.append(z_coil2)
    else:
        for cell_index in range(n_cells):
            z_0 = (-n_cells/2 + cell_index) * cell_length
            z_coil1 = z_0 - cell_length/2 + coil_cell_z
            z_coil2 = z_0 + cell_length/2 - coil_cell_z
            coil_global_z.append(z_coil1)
            coil_global_z.append(z_coil2)
    return coil_global_z


def build_periodic_channel(N_cells, L_cell, z_coils, polarities=None):
    """Build coil positions for a periodic channel centered at z=0.

    Each cell contains coils at z_coils[i] and their mirrors at (L_cell - z_coils[i]).
    Mirror coils get flipped polarity.

    Parameters
    ----------
    N_cells    : int         - Number of cells
    L_cell     : float       - Cell length [m]
    z_coils    : list[float] - Coil positions within cell, measured from cell start [m]
                               (only the "left" coils; mirrors are generated automatically)
    polarities : list[int], optional - Polarity (+1/-1) for each coil in z_coils.
                               Mirrors get opposite sign. Default: all +1.

    Returns
    -------
    list of dicts: [{'z': global_z, 'polarity': +1/-1, 'coil_type': index}, ...]
                   sorted by z position
    """
    if polarities is None:
        polarities = [1] * len(z_coils)

    z_channel_start = -N_cells * L_cell / 2

    coils = []
    for cell_idx in range(N_cells):
        z_cell_start = z_channel_start + cell_idx * L_cell

        for coil_idx, (z_local, pol) in enumerate(zip(z_coils, polarities)):
            coils.append({
                'z': z_cell_start + z_local,
                'polarity': pol,
                'coil_type': coil_idx
            })
            coils.append({
                'z': z_cell_start + L_cell - z_local,
                'polarity': -pol,
                'coil_type': coil_idx
            })

    coils.sort(key=lambda c: c['z'])
    return coils


def build_coil_beamline(N_cells=None, L_cell=None, z_coils=None, coil_templates=None,
                        polarities=None, fixed=True, name="Channel",
                        channel: MuonCoolingChannel = None):
    """Build a MuonCoolingChannel from a periodic channel definition.

    Parameters
    ----------
    N_cells        : int         - Number of cells
    L_cell         : float       - Cell length [m]
    z_coils        : list[float] - Coil positions within cell from cell start [m]
    coil_templates : list[dict]  - One per coil type
    polarities     : list[int]   - Polarity per coil type
    fixed          : bool        - Whether channel coils are fixed for matching
    name           : str         - Beamline name (if channel is None)
    channel        : MuonCoolingChannel - Existing channel to add elements to

    Returns
    -------
    MuonCoolingChannel object with all channel coils
    """
    if channel is not None:
        if L_cell is not None and not np.isclose(L_cell, channel.cell_length):
            raise ValueError(f"L_cell {L_cell} != channel.cell_length {channel.cell_length}")
        L_cell = channel.cell_length
        if N_cells is not None and N_cells > channel.n_cells:
            raise ValueError(f"N_cells {N_cells} != channel.n_cells {channel.n_cells}")
        N_cells = channel.n_cells
    elif L_cell is None or N_cells is None:
        raise ValueError("Must provide N_cells and L_cell if channel is None")

    positions = build_periodic_channel(N_cells, L_cell, z_coils, polarities)

    if channel is None:
        channel = MuonCoolingChannel(n_cells=N_cells, cell_length=L_cell, total_width=1.0, name=name)

    for i, pos in enumerate(positions):
        tmpl = coil_templates[pos['coil_type']]
        pol = pos['polarity']

        J = tmpl['currDensity'] * pol

        coil = SolenoidCoil(
            z_center=pos['z'],
            r_in=tmpl['r_in'],
            r_thick=tmpl['r_thick'],
            N_pancakes=tmpl['N_pancakes'],
            L_pancake=tmpl['L_pancake'],
            L_spacing=tmpl['L_spacing'],
            currDensity=J,
            N_sheets=tmpl.get('N_sheets', 10),
            fixed=fixed,
            name=f"{tmpl['name']}_{i}",
            material=tmpl.get('material', 'G4_Cu'),
        )
        channel.add_element(coil)

    return channel


def build_dipole_beamline(n_cells=None, cell_length=None, coil_cell_z=None, dipole_template=None,
                          polarities=None, fixed=True, name="DipoleChannel",
                          channel: MuonCoolingChannel = None):
    """Build a MuonCoolingChannel of Dipoles using coil_placement_z for positions.

    Parameters
    ----------
    n_cells         : int      - Number of cells
    cell_length     : float    - Cell length [m]
    coil_cell_z     : float    - Dipole offset from cell edge [m]
    dipole_template : dict     - Template for dipole properties
    polarities      : list[int] - Polarity pattern
    fixed           : bool     - Whether dipoles are fixed
    name            : str      - Beamline name (if channel is None)
    channel         : MuonCoolingChannel - Existing channel to add elements to

    Returns
    -------
    MuonCoolingChannel object with dipoles added
    """
    if channel is not None:
        if cell_length is not None and not np.isclose(cell_length, channel.cell_length):
            raise ValueError(f"cell_length {cell_length} != channel.cell_length {channel.cell_length}")
        cell_length = channel.cell_length
        if n_cells is not None and n_cells > channel.n_cells:
            raise ValueError(f"n_cells {n_cells} != channel.n_cells {channel.n_cells}")
        n_cells = channel.n_cells
    elif cell_length is None or n_cells is None:
        raise ValueError("Must provide n_cells and cell_length if channel is None")

    if channel is None:
        channel = MuonCoolingChannel(n_cells=n_cells, cell_length=cell_length, total_width=1.0, name=name,
                                     polarities=polarities)

    positions = coil_placement_z(n_cells, cell_length, coil_cell_z)

    for i, z in enumerate(positions):
        pol = channel.polarities[i % len(channel.polarities)]
        dipole = Dipole(
            z_center=z,
            field_strength=dipole_template['field_strength'] * pol,
            aperture=dipole_template.get('aperture', 0.2),
            length_z=dipole_template.get('length_z', 0.1),
            enge_coefficient=dipole_template.get('enge_coefficient', 5.5),
            fixed=fixed,
            name=f"{dipole_template.get('name', 'Dipole')}_{i}",
        )
        channel.add_element(dipole)

    return channel


def build_absorber_beamline(n_cells=None, cell_length=None, absorber_template=None,
                            fixed=True, name="AbsorberChannel", channel: MuonCoolingChannel = None):
    """Place one absorber per cell in a MuonCoolingChannel.

    Wedge rotation is derived from channel.polarities: positive polarity → math.pi (absdown),
    negative → 0.0 (absup). Polarities are sampled with a stride of
    len(polarities) // n_cells so that one absorber covers each polarity pair.

    Parameters
    ----------
    n_cells           : int         - Number of cells
    cell_length       : float       - Cell length [m]
    absorber_template : dict        - Template for absorber properties
    fixed             : bool        - Whether absorbers are fixed
    name              : str         - Beamline name (if channel is None)
    channel           : MuonCoolingChannel - Existing channel to add elements to
    """
    if channel is not None:
        if cell_length is not None and not np.isclose(cell_length, channel.cell_length):
            raise ValueError(f"cell_length {cell_length} != channel.cell_length {channel.cell_length}")
        cell_length = channel.cell_length
        if n_cells is not None and n_cells > channel.n_cells:
            raise ValueError(f"n_cells {n_cells} != channel.n_cells {channel.n_cells}")
        n_cells = channel.n_cells
    elif cell_length is None or n_cells is None:
        raise ValueError("Must provide n_cells and cell_length if channel is None")

    if channel is None:
        channel = MuonCoolingChannel(n_cells=n_cells, cell_length=cell_length, total_width=1.0, name=name)

    stride  = max(1, len(channel.polarities) // n_cells)
    z_start = -n_cells * cell_length / 2

    for cell_idx in range(n_cells):
        cell_z = z_start + cell_idx * cell_length
        pol    = channel.polarities[(cell_idx * stride) % len(channel.polarities)]
        angle  = math.pi if pol > 0 else 0.0

        channel.add_element(Absorber(
            z_center=cell_z,
            wedge_rotation_angle=angle,
            **{k: v for k, v in absorber_template.items()
               if k not in ('wedge_rotation_angle', 'name')},
            fixed=fixed,
            name=f"{absorber_template.get('name', 'Absorber')}_{cell_idx}",
        ))

    return channel


def build_rf_beamline(n_cells=None, cell_length=None, n_rf_cells=None, rf_spacing=None,
                      rf_template=None, fixed=True, name="RFChannel", channel: MuonCoolingChannel = None):
    """Place n_rf_cells RF cavities per cell in a MuonCoolingChannel.

    Parameters
    ----------
    n_cells     : int      - Number of cells
    cell_length : float    - Cell length [m]
    n_rf_cells  : int      - Number of RF cavities per cell
    rf_spacing  : float    - Spacing between RF cavities [m]
    rf_template : dict     - Template for RF properties
    fixed       : bool     - Whether cavities are fixed
    name        : str      - Beamline name (if channel is None)
    channel     : MuonCoolingChannel - Existing channel to add elements to
    """
    if channel is not None:
        if cell_length is not None and not np.isclose(cell_length, channel.cell_length):
            raise ValueError(f"cell_length {cell_length} != channel.cell_length {channel.cell_length}")
        cell_length = channel.cell_length
        if n_cells is not None and n_cells > channel.n_cells:
            raise ValueError(f"n_cells {n_cells} != channel.n_cells {channel.n_cells}")
        n_cells = channel.n_cells
    elif cell_length is None or n_cells is None:
        raise ValueError("Must provide n_cells and cell_length if channel is None")

    if channel is None:
        channel = MuonCoolingChannel(n_cells=n_cells, cell_length=cell_length, total_width=1.0, name=name)

    z_start  = -n_cells * cell_length / 2
    rf_count = 0

    for cell_idx in range(n_cells):
        cell_z    = z_start + cell_idx * cell_length
        rf_center = cell_z + cell_length / 2

        for rf_idx in range(1, n_rf_cells + 1):
            z_rf = rf_center + (rf_idx - 0.5 - n_rf_cells / 2) * rf_spacing
            channel.add_element(RFCavity(
                z_center=z_rf,
                **{k: v for k, v in rf_template.items() if k != 'name'},
                fixed=fixed,
                name=f"{rf_template.get('name', 'RF')}_{rf_count}",
            ))
            rf_count += 1

    return channel
