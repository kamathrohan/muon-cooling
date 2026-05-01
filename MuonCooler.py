import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional
import pandas as pd

class BeamlineElement:
    """Base class for beamline elements."""
    
    def __init__(self, z_center: float, name: Optional[str] = None):
        self.z_center = z_center
        self.name = name or f"Element@{z_center:.3f}m"
    
    def get_field(self, z: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Subclasses must implement get_field()")
    
    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}', z_center={self.z_center})"


class SolenoidCoil(BeamlineElement):
    """Solenoid coil modeled as stacked cylindrical current sheets.
    
    The coil can be a single block or an array of pancakes with spacings.
    Current can be specified via current density (preferred) or total current.
    
    Parameters
    ----------
    z_center     : float   - Axial center position [m]
    r_in         : float   - Inner radius [m]
    r_thick      : float   - Radial thickness [m]
    N_sheets     : int     - Number of radial current sheets per pancake
    fixed        : bool    - If True, current is held fixed during matching
    name         : str     - Name for this coil
    
    Current (specify one):
        currDensity : float - Current density [A/mm^2] (preferred)
        Itot        : float - Total current [A] (alternative)
    
    Geometry (two modes):
        Simple block:
            L : float - Axial length [m]
        Pancake array:
            N_pancakes : int   - Number of pancakes
            L_pancake  : float - Axial length of each pancake [m]
            L_spacing  : float - Gap between adjacent pancakes [m]
    """
    
    def __init__(self, z_center: float, r_in: float, r_thick: float,
                 L: float = None,
                 N_pancakes: int = 1, L_pancake: float = None, L_spacing: float = 0.0,
                 currDensity: float = None, Itot: float = None,
                 N_sheets: int = 10, fixed: bool = False, name: Optional[str] = None,
                 material: str = "G4_Cu"):
        super().__init__(z_center, name)
        self.r_in = r_in
        self.r_thick = r_thick
        self.N_sheets = N_sheets
        self.fixed = fixed
        self.material = material
        
        # --- Geometry ---
        if L is not None and L_pancake is not None:
            raise ValueError("Specify either L (single block) or L_pancake/N_pancakes (pancake array), not both.")
        
        if L is not None:
            # Simple block mode: treat as 1 pancake
            self.N_pancakes = 1
            self.L_pancake = L
            self.L_spacing = 0.0
        elif L_pancake is not None:
            self.N_pancakes = N_pancakes
            self.L_pancake = L_pancake
            self.L_spacing = L_spacing
        else:
            raise ValueError("Must specify either L or L_pancake.")
        
        # Total coil length
        self.L_total = self.N_pancakes * self.L_pancake + (self.N_pancakes - 1) * self.L_spacing
        
        # --- Current ---
        if currDensity is not None and Itot is not None:
            raise ValueError("Specify either currDensity or Itot, not both.")
        
        if currDensity is not None:
            self.currDensity = currDensity  # A/mm^2
            # Itot per pancake = J * L_pancake[mm] * r_thick[mm]
            self._Itot_per_pancake = currDensity * (self.L_pancake * 1e3) * (self.r_thick * 1e3)
        elif Itot is not None:
            self._Itot_per_pancake = Itot / self.N_pancakes
            self.currDensity = self._Itot_per_pancake / (self.L_pancake * 1e3) / (self.r_thick * 1e3)
        else:
            raise ValueError("Must specify either currDensity or Itot.")
    
    @property
    def Itot(self):
        """Total current across all pancakes [A]."""
        return self._Itot_per_pancake * self.N_pancakes
    
    @Itot.setter
    def Itot(self, value):
        """Set total current; redistributes equally across pancakes."""
        self._Itot_per_pancake = value / self.N_pancakes
        self.currDensity = self._Itot_per_pancake / (self.L_pancake * 1e3) / (self.r_thick * 1e3)
    
    @property
    def currDensity_value(self):
        """Current density [A/mm^2]."""
        return self.currDensity
    
    @currDensity_value.setter 
    def currDensity_value(self, value):
        """Set current density; updates all pancake currents."""
        self.currDensity = value
        self._Itot_per_pancake = value * (self.L_pancake * 1e3) * (self.r_thick * 1e3)
        
    def _bzSheet(self, z: np.ndarray, Itot: float, L: float, 
                 r: float, z0: float) -> np.ndarray:
        """Field from a single cylindrical current sheet."""
        plus = z - z0 + L / 2
        minus = z - z0 - L / 2
        mu0 = 4 * np.pi * 1e-7
        return (mu0 * Itot / (2 * L)) * (
            plus / np.sqrt(plus**2 + r**2) - 
            minus / np.sqrt(minus**2 + r**2)
        )
    
    def _pancake_centers(self) -> np.ndarray:
        """Compute the z-center of each pancake."""
        pitch = self.L_pancake + self.L_spacing
        # Center the pancake array on self.z_center
        start = self.z_center - (self.N_pancakes - 1) * pitch / 2
        return np.array([start + i * pitch for i in range(self.N_pancakes)])
    
    def get_field(self, z: np.ndarray) -> np.ndarray:
        bz = np.zeros_like(z, dtype=float)
        dr = self.r_thick / self.N_sheets
        Itot_per_sheet = self._Itot_per_pancake / self.N_sheets
        
        for zp in self._pancake_centers():
            for i in range(self.N_sheets):
                r_sheet = self.r_in + (i + 0.5) * dr
                bz += self._bzSheet(z, Itot_per_sheet, self.L_pancake, r_sheet, zp)
        return bz
    
    def __repr__(self):
        polarity = "+" if self._Itot_per_pancake >= 0 else "-"
        tag = " [FIXED]" if self.fixed else ""
        if self.N_pancakes == 1:
            geom = f"L={self.L_pancake}m"
        else:
            geom = f"{self.N_pancakes}x{self.L_pancake*1e3:.1f}mm (gap={self.L_spacing*1e3:.1f}mm)"
        return (f"SolenoidCoil('{self.name}', z={self.z_center:.3f}m, "
                f"J={polarity}{abs(self.currDensity):.1f}A/mm², {geom}){tag}")


class Dipole(BeamlineElement):
    def __init__(self, z_center: float, field_strength: float,
                 aperture: float = 0.2, length_z: float = 0.1,
                 enge_coefficient: float = 5.5,
                 fixed: bool = False, name: Optional[str] = None):
        super().__init__(z_center, name)
        self.field_strength   = field_strength
        self.aperture         = aperture
        self.length_z         = length_z
        self.enge_coefficient = enge_coefficient
        self.fixed            = fixed

    def get_field(self, z: np.ndarray) -> np.ndarray:
        return np.zeros_like(z, dtype=float)

    def __repr__(self):
        tag = " [FIXED]" if self.fixed else ""
        return (f"Dipole('{self.name}', z={self.z_center:.3f}m, "
                f"B={self.field_strength:.3f}T){tag}")


class Absorber(BeamlineElement):
    def __init__(self, z_center: float, absorber_type: str = "wedge",
                 material: str = "G4_LITHIUM_HYDRIDE",
                 cylinder_length: float = 1.10, cylinder_radius: float = 0.24,
                 wedge_opening_angle: float = 0.17453292519943295,
                 wedge_height: float = 0.24,
                 wedge_rotation_angle: float = 0.0,
                 wedge_offset_x: float = 0.0, wedge_offset_y: float = 0.0,
                 wedge_apex_to_base: float = 0.1134256,
                 fixed: bool = False, name: Optional[str] = None):
        super().__init__(z_center, name)
        self.absorber_type        = absorber_type
        self.material             = material
        self.cylinder_length      = cylinder_length
        self.cylinder_radius      = cylinder_radius
        self.wedge_opening_angle  = wedge_opening_angle
        self.wedge_height         = wedge_height
        self.wedge_rotation_angle = wedge_rotation_angle
        self.wedge_offset_x       = wedge_offset_x
        self.wedge_offset_y       = wedge_offset_y
        self.wedge_apex_to_base   = wedge_apex_to_base
        self.fixed                = fixed

    def get_field(self, z: np.ndarray) -> np.ndarray:
        return np.zeros_like(z, dtype=float)

    def __repr__(self):
        tag = " [FIXED]" if self.fixed else ""
        return (f"Absorber('{self.name}', z={self.z_center:.3f}m, "
                f"type={self.absorber_type}, mat={self.material}){tag}")


class RFCavity(BeamlineElement):
    def __init__(self, z_center: float, time_offset: float = 0.0,
                 length: float = 0.18856, voltage: float = 30.0,
                 phase: float = -1.5707963267948966 + 0.3490658503988659,
                 frequency: float = 704e6,
                 window_thickness: float = 0.0, window_material: str = "G4_Be",
                 window_radius: float = 0.0816,
                 cavity_material: str = "G4_Cu",
                 cavity_vacuum_material: str = "vacuum",
                 cavity_radius: float = 0.163, cavity_thickness: float = 0.003,
                 fixed: bool = False, name: Optional[str] = None):
        super().__init__(z_center, name)
        self.time_offset            = time_offset
        self.length                 = length
        self.voltage                = voltage
        self.phase                  = phase
        self.frequency              = frequency
        self.window_thickness       = window_thickness
        self.window_material        = window_material
        self.window_radius          = window_radius
        self.cavity_material        = cavity_material
        self.cavity_vacuum_material = cavity_vacuum_material
        self.cavity_radius          = cavity_radius
        self.cavity_thickness       = cavity_thickness
        self.fixed                  = fixed

    def get_field(self, z: np.ndarray) -> np.ndarray:
        return np.zeros_like(z, dtype=float)

    def __repr__(self):
        tag = " [FIXED]" if self.fixed else ""
        return (f"RFCavity('{self.name}', z={self.z_center:.3f}m, "
                f"V={self.voltage}MV, f={self.frequency/1e6:.0f}MHz){tag}")


class Beamline:
    """Container for beamline elements with field calculation."""
    
    def __init__(self, name: str = "Beamline"):
        self.name = name
        self.elements: List[BeamlineElement] = []
    
    def add_element(self, element: BeamlineElement) -> None:
        self.elements.append(element)
        
    def add_elements(self, elements: List[BeamlineElement]) -> None:
        self.elements.extend(elements)
    
    def remove_element(self, index: int) -> BeamlineElement:
        return self.elements.pop(index)
    
    def clear(self) -> None:
        self.elements.clear()
    
    def get_field(self, z: np.ndarray) -> np.ndarray:
        bz_total = np.zeros_like(z, dtype=float)
        for element in self.elements:
            bz_total += element.get_field(z)
        return bz_total
    
    def get_element_fields(self, z: np.ndarray) -> Dict[str, np.ndarray]:
        fields = {}
        for element in self.elements:
            fields[element.name] = element.get_field(z)
        return fields
    
    def plot_field(self, z_range: tuple = None, n_points: int = 1000,
                   show_elements: bool = True, **plot_kwargs):
        if z_range is None:
            if not self.elements:
                raise ValueError("No elements in beamline")
            z_positions = [elem.z_center for elem in self.elements]
            margin = 1.0
            z_range = (min(z_positions) - margin, max(z_positions) + margin)
        
        z = np.linspace(z_range[0], z_range[1], n_points)
        bz = self.get_field(z)
        
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(z, bz, **plot_kwargs)
        
        if show_elements:
            for elem in self.elements:
                if isinstance(elem, SolenoidCoil):
                    color = 'blue' if elem._Itot_per_pancake > 0 else 'red'
                    ax.axvline(elem.z_center, color=color, linestyle='--', 
                              alpha=0.3, linewidth=1)
        
        ax.set_xlabel("z [m]", fontsize=14)
        ax.set_ylabel(r"$B_z$ [T]", fontsize=14)
        ax.grid(True)
        plt.tight_layout()
        return fig, ax
    
    def __len__(self):
        return len(self.elements)
    
    def __getitem__(self, index):
        return self.elements[index]
    
    def __repr__(self):
        return f"Beamline(name='{self.name}', n_elements={len(self.elements)})"

    def compute_rf_time_offsets(self, momentum_MeV_c: float,
                                mass_MeV_c: float = 105.66,
                                z_ref: float = None) -> Dict[str, float]:
        """Compute and set time_offset on every RFCavity from transit time.

        t_i = (z_i - z_ref) / (beta * c)

        Parameters
        ----------
        momentum_MeV_c : float - Reference particle momentum [MeV/c]
        mass_MeV_c     : float - Particle rest mass [MeV/c²] (default: muon)
        z_ref          : float - Reference z [m]; defaults to leftmost element.

        Returns
        -------
        dict mapping cavity name -> time_offset [ns]
        """
        import math
        C_M_PER_NS = 0.299792458  # speed of light [m/ns]

        E    = math.sqrt(momentum_MeV_c**2 + mass_MeV_c**2)
        beta = momentum_MeV_c / E

        if z_ref is None:
            z_ref = min(e.z_center for e in self.elements)

        offsets = {}
        for elem in self.elements:
            if isinstance(elem, RFCavity):
                elem.time_offset = (elem.z_center - z_ref) / (beta * C_M_PER_NS)
                offsets[elem.name] = elem.time_offset

        return offsets

    def summary(self):
        print(f"\n{self.name}")
        print("=" * 60)
        print(f"Total elements: {len(self.elements)}\n")
        for i, elem in enumerate(self.elements):
            print(f"[{i}] {elem}")

    def build_pancake_dataframe(self) -> pd.DataFrame:
        """Return a DataFrame with one row per sub-pancake across all SolenoidCoils."""
        rows = []
        for coil_idx, elem in enumerate(self.elements):
            if not isinstance(elem, SolenoidCoil):
                continue
            centers = elem._pancake_centers()
            polarity = 1 if elem._Itot_per_pancake >= 0 else -1
            for p_idx, z_p in enumerate(centers):
                rows.append({
                    'coil_name':     f"{elem.name}_p{p_idx}",
                    'coil_index':    coil_idx,
                    'pancake_index': p_idx,
                    'z_center':      z_p,
                    'z_coil_center': elem.z_center,
                    'Itot_pancake':  elem._Itot_per_pancake,
                    'currDensity':   elem.currDensity,
                    'polarity':      polarity,
                    'r_in':          elem.r_in,
                    'r_thick':       elem.r_thick,
                    'L_pancake':     elem.L_pancake,
                    'fixed':         elem.fixed,
                    'material':      elem.material,
                    'coil_offset_x': 0.0,
                    'coil_offset_y': 0.0,
                    'coil_tilt_x':   0.0,
                    'coil_tilt_y':   0.0,
                    'coil_tilt_z':   0.0,
                })
        return pd.DataFrame(rows)


class MuonCoolingChannel(Beamline):
    def __init__(
        self,
        n_cells: int,
        cell_length: float,
        total_length: float,
        total_width: float,
        name: str = "MuonCoolingChannel",
        magnetic_field_model: str = "solenoidsheet",
        magnetic_field_method: str = "grid",
        dipole_field_model: str = "dipole",
        interpolator: str = "linear",
        electric_field_model: str = "rfpillbox",
        reference_momentum: float = 200.0,
        on_axis_tolerance: float = 2e-2,
    ):
        super().__init__(name)
        self.n_cells = n_cells
        self.cell_length = cell_length
        self.total_length = total_length
        self.total_width = total_width

        # Rendering / BDSIM field settings
        self.magnetic_field_model = magnetic_field_model
        self.magnetic_field_method = magnetic_field_method
        self.dipole_field_model = dipole_field_model
        self.interpolator = interpolator
        self.electric_field_model = electric_field_model
        self.reference_momentum = reference_momentum
        self.on_axis_tolerance = on_axis_tolerance

    def __repr__(self):
        return (f"MuonCoolingChannel(name='{self.name}', n_cells={self.n_cells}, "
                f"cell_length={self.cell_length}m, total_length={self.total_length}m, "
                f"total_width={self.total_width}m, n_elements={len(self.elements)})")

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

print("Beamline classes loaded.")


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


def build_channel_beamline(N_cells=None, L_cell=None, z_coils=None, coil_templates=None, 
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
        # Default width if creating new channel
        channel = MuonCoolingChannel(n_cells=N_cells, cell_length=L_cell, total_width=1.0, name=name)
    
    for i, pos in enumerate(positions):
        tmpl = coil_templates[pos['coil_type']]
        pol = pos['polarity']
        
        # Compute Itot from template's currDensity, apply polarity
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

    if polarities is None:
        polarities = [1, -1, -1, 1]

    if channel is None:
        channel = MuonCoolingChannel(n_cells=n_cells, cell_length=cell_length, total_width=1.0, name=name)

    positions = coil_placement_z(n_cells, cell_length, coil_cell_z)

    for i, z in enumerate(positions):
        pol = polarities[i % len(polarities)]
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
                             rotation_angles=None, fixed=True,
                             name="AbsorberChannel", channel: MuonCoolingChannel = None):
    """Place one absorber per cell in a MuonCoolingChannel.

    Parameters
    ----------
    n_cells           : int         - Number of cells
    cell_length       : float       - Cell length [m]
    absorber_template : dict        - Template for absorber properties
    rotation_angles   : list[float] - Wedge rotation angles cycled per cell
    fixed             : bool        - Whether absorbers are fixed
    name              : str         - Beamline name (if channel is None)
    channel           : MuonCoolingChannel - Existing channel to add elements to
    """
    import math
    if channel is not None:
        if cell_length is not None and not np.isclose(cell_length, channel.cell_length):
            raise ValueError(f"cell_length {cell_length} != channel.cell_length {channel.cell_length}")
        cell_length = channel.cell_length
        if n_cells is not None and n_cells > channel.n_cells:
            raise ValueError(f"n_cells {n_cells} != channel.n_cells {channel.n_cells}")
        n_cells = channel.n_cells
    elif cell_length is None or n_cells is None:
        raise ValueError("Must provide n_cells and cell_length if channel is None")

    if rotation_angles is None:
        rotation_angles = [math.pi, 0.0]

    if channel is None:
        channel = MuonCoolingChannel(n_cells=n_cells, cell_length=cell_length, total_width=1.0, name=name)

    z_start = -n_cells * cell_length / 2

    for cell_idx in range(n_cells):
        cell_z = z_start + cell_idx * cell_length
        angle  = rotation_angles[cell_idx % len(rotation_angles)]

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





def generate_sampler_lines(n_samplers: int,
                           total_length_m: float,
                           aper_m: float = 5.0,
                           reference_element: str = "mc1",
                           reference_element_number: int = 0) -> str:
    """Generate BDSIM samplerplacement lines evenly spaced across the channel.

    Parameters
    ----------
    n_samplers              : int   - Number of samplers to place
    total_length_m          : float - Full channel length [m]; samplers span ±half
    aper_m                  : float - Rectangular aperture half-size [m]
    reference_element       : str   - BDSIM element name to reference
    reference_element_number: int   - referenceElementNumber value

    Returns
    -------
    str - Multi-line BDSIM sampler placement block
    """
    half_mm = total_length_m / 2 * 1e3
    positions = np.linspace(-half_mm, half_mm, n_samplers)
    lines = []
    for i, s in enumerate(positions, 1):
        lines.append(
            f's{i}: samplerplacement, referenceElement="{reference_element}", '
            f'referenceElementNumber={reference_element_number}, '
            f's={s:.1f}*mm, shape="rectangular",'
            f'aper1={aper_m}*m, aper2={aper_m}*m;'
        )
    return '\n'.join(lines)


def generate_beam_block(mode: str,
                        particle: str = "mu+",
                        momentum_MeV: float = 200.0,
                        x0_mm: float = 0.0,
                        y0_mm: float = 0.0,
                        z0_m: float = 0.0,
                        distr_type: str = "userfile",
                        distr_file: str = "",
                        distr_file_format: str = "t[ns]:x[m]:y[m]:z[m]:xp[rad]:yp[rad]:zp[rad]:-:E[GeV]",
                        nlines_ignore: int = 1) -> str:
    """Generate a BDSIM beam block.

    Parameters
    ----------
    mode             : str   - One of 'reference', 'offset', or 'beam'
                               'reference' : particle + momentum only
                               'offset'    : particle + momentum + X0/Y0/Z0
                               'beam'      : particle + momentum + distribution file
    particle         : str   - Particle type (default "mu+")
    momentum_MeV     : float - Reference momentum [MeV/c]
    x0_mm, y0_mm     : float - Transverse offsets [mm] (offset mode)
    z0_m             : float - Longitudinal offset [m] (offset mode)
    distr_type       : str   - BDSIM distribution type (beam mode)
    distr_file       : str   - Path to distribution file (beam mode)
    distr_file_format: str   - Column format string (beam mode)
    nlines_ignore    : int   - Header lines to skip in distribution file (beam mode)

    Returns
    -------
    str - BDSIM beam command block
    """
    if mode == "reference":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV;'
        )
    elif mode == "offset":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV,\n'
            f'      X0={x0_mm}*mm,\n'
            f'      Y0={y0_mm}*mm,\n'
            f'      Z0={z0_m}*m;'
        )
    elif mode == "beam":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV,\n'
            f'      distrType="{distr_type}",\n'
            f'      distrFile="{distr_file}",\n'
            f'      distrFileFormat = "{distr_file_format}",\n'
            f'      nlinesIgnore={nlines_ignore};'
        )
    else:
        raise ValueError(f"Unknown beam mode '{mode}'. Choose 'reference', 'offset', or 'beam'.")


MUON_MASS_MEV_C2 = 105.66

def render_gmad(beamline: MuonCoolingChannel,
                tpl_path: str,
                out_gmad: str,
                n_samplers: int = 129,
                sampler_aper_m: float = 5.0,
                beam_mode: str = "beam",
                beam_kwargs: dict = None) -> None:
    """Render a full channel gmad from channel.tpl and a MuonCoolingChannel object.

    Parameters
    ----------
    beamline : MuonCoolingChannel - channel with Solenoid/Dipole/Absorber/RFCavity elements
    tpl_path : str          - path to channel.tpl
    out_gmad : str          - output gmad path
    """
    from jinja2 import Template

    beamline.compute_rf_time_offsets(beamline.reference_momentum, MUON_MASS_MEV_C2)

    pancake_df = beamline.build_pancake_dataframe()
    df = pancake_df.sort_values(['coil_index', 'pancake_index']).reset_index(drop=True)

    dipoles   = [e for e in beamline.elements if isinstance(e, Dipole)]
    absorbers = [e for e in beamline.elements if isinstance(e, Absorber)]
    rf_cavs   = [e for e in beamline.elements if isinstance(e, RFCavity)]

    def fmt_array(values, fmt='.6g'):
        formatted = [format(v, fmt) for v in values]
        if len(set(formatted)) == 1:
            return '{' + formatted[0] + '}'
        return '{' + ', '.join(formatted) + '}'

    def fmt_mm_array(values_m):
        formatted = [f'{v*1e3:.4f}*mm' for v in values_m]
        if len(set(formatted)) == 1:
            return '{' + formatted[0] + '}'
        return '{' + ', '.join(formatted) + '}'

    def fmt_str_array(values):
        if len(set(values)) == 1:
            return '{"' + values[0] + '"}'
        return '{' + ', '.join(f'"{v}"' for v in values) + '}'

    n_coils = len(df)
    context = {
        # channel properties
        'total_length':      beamline.total_length,
        'total_width':       beamline.total_width,
        'on_axis_tolerance': beamline.on_axis_tolerance,
        'coil_material':     fmt_str_array(list(df['material'])),

        # channel field models
        'magnetic_field_model':  beamline.magnetic_field_model,
        'magnetic_field_method': beamline.magnetic_field_method,
        'dipole_field_model':    beamline.dipole_field_model,
        'interpolator':          beamline.interpolator,
        'electric_field_model':  beamline.electric_field_model,

        # coils
        'n_coils':               n_coils,
        'coil_inner_radius':     fmt_mm_array(df['r_in']),
        'coil_radial_thickness': fmt_mm_array(df['r_thick']),
        'coil_length_z':         fmt_array(df['L_pancake']),
        'coil_current':          fmt_array(df['currDensity']),
        'coil_offset_x':         fmt_array(df['coil_offset_x']),
        'coil_offset_y':         fmt_array(df['coil_offset_y']),
        'coil_offset_z':         fmt_array(df['z_center'], '.9g'),
        'coil_tilt_x':           fmt_array(df['coil_tilt_x']),
        'coil_tilt_y':           fmt_array(df['coil_tilt_y']),
        'coil_tilt_z':           fmt_array(df['coil_tilt_z']),

        # dipoles
        'n_dipoles':              len(dipoles),
        'dipole_aperture':        fmt_mm_array([d.aperture      for d in dipoles]),
        'dipole_length_z':        fmt_mm_array([d.length_z      for d in dipoles]),
        'dipole_field_strength':  fmt_array([d.field_strength   for d in dipoles]),
        'dipole_enge_coefficient':fmt_array([d.enge_coefficient for d in dipoles]),
        'dipole_offset_z':        fmt_array([d.z_center         for d in dipoles], '.9g'),

        # absorbers
        'n_absorbers':                len(absorbers),
        'absorber_type':              fmt_str_array([a.absorber_type       for a in absorbers]),
        'absorber_material':          fmt_str_array([a.material            for a in absorbers]),
        'absorber_offset_z':          fmt_array([a.z_center               for a in absorbers], '.9g'),
        'absorber_cylinder_length':   fmt_mm_array([a.cylinder_length      for a in absorbers]),
        'absorber_cylinder_radius':   fmt_mm_array([a.cylinder_radius      for a in absorbers]),
        'absorber_wedge_opening_angle': fmt_array([a.wedge_opening_angle   for a in absorbers]),
        'absorber_wedge_height':      fmt_mm_array([a.wedge_height         for a in absorbers]),
        'absorber_wedge_rotation_angle': fmt_array([a.wedge_rotation_angle for a in absorbers]),
        'absorber_wedge_offset_x':    fmt_array([a.wedge_offset_x          for a in absorbers]),
        'absorber_wedge_offset_y':    fmt_array([a.wedge_offset_y          for a in absorbers]),
        'absorber_wedge_apex_to_base':fmt_mm_array([a.wedge_apex_to_base   for a in absorbers]),

        # rf cavities
        'n_rf_cavities':          len(rf_cavs),
        'rf_offset_z':            fmt_array([r.z_center          for r in rf_cavs], '.9g'),
        'rf_time_offset':         fmt_array([r.time_offset        for r in rf_cavs]),
        'rf_length':              fmt_mm_array([r.length          for r in rf_cavs]),
        'rf_voltage':             fmt_array([r.voltage            for r in rf_cavs]),
        'rf_phase':               fmt_array([r.phase              for r in rf_cavs]),
        'rf_frequency':           fmt_array([r.frequency          for r in rf_cavs]),
        'rf_window_thickness':    fmt_mm_array([r.window_thickness for r in rf_cavs]),
        'rf_window_material':     fmt_str_array([r.window_material for r in rf_cavs]),
        'rf_window_radius':       fmt_mm_array([r.window_radius   for r in rf_cavs]),
        'rf_cavity_material':     fmt_str_array([r.cavity_material for r in rf_cavs]),
        'rf_cavity_vacuum_material': fmt_str_array([r.cavity_vacuum_material for r in rf_cavs]),
        'rf_cavity_radius':       fmt_mm_array([r.cavity_radius   for r in rf_cavs]),
        'rf_cavity_thickness':    fmt_mm_array([r.cavity_thickness for r in rf_cavs]),
    }

    context['sampler_lines'] = generate_sampler_lines(
        n_samplers, beamline.total_length, aper_m=sampler_aper_m
    )
    context['beam_block'] = generate_beam_block(
        beam_mode, momentum_MeV=beamline.reference_momentum, **(beam_kwargs or {})
    )

    with open(tpl_path) as f:
        rendered = Template(f.read()).render(**context)

    with open(out_gmad, 'w') as f:
        f.write(rendered)

    print(f"Written {out_gmad}: {n_coils} coils, {len(dipoles)} dipoles, "
          f"{len(absorbers)} absorbers, {len(rf_cavs)} RF cavities, "
          f"{n_samplers} samplers")



#TODO: Error handling of muon cooling channel size < z positions
