import numpy as np
from typing import Optional


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
