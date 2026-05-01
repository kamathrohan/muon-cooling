import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Optional

from .elements import BeamlineElement, SolenoidCoil, RFCavity, Absorber, Dipole


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

    def build_rf_dataframe(self, global_z_offset: float = 0.0) -> pd.DataFrame:
        """Return a DataFrame with one row per RFCavity."""
        rows = []
        for elem in self.elements:
            if not isinstance(elem, RFCavity):
                continue
            rows.append({
                'name':                   elem.name,
                'z_local':                elem.z_center,
                'z_global':               elem.z_center + global_z_offset,
                'time_offset':            elem.time_offset,
                'length':                 elem.length,
                'voltage':                elem.voltage,
                'phase':                  elem.phase,
                'frequency':              elem.frequency,
                'window_thickness':       elem.window_thickness,
                'window_material':        elem.window_material,
                'window_radius':          elem.window_radius,
                'cavity_material':        elem.cavity_material,
                'cavity_vacuum_material': elem.cavity_vacuum_material,
                'cavity_radius':          elem.cavity_radius,
                'cavity_thickness':       elem.cavity_thickness,
                'fixed':                  elem.fixed,
            })
        return pd.DataFrame(rows)

    def build_absorber_dataframe(self, global_z_offset: float = 0.0) -> pd.DataFrame:
        """Return a DataFrame with one row per Absorber."""
        rows = []
        for elem in self.elements:
            if not isinstance(elem, Absorber):
                continue
            rows.append({
                'name':                 elem.name,
                'z_local':              elem.z_center,
                'z_global':             elem.z_center + global_z_offset,
                'absorber_type':        elem.absorber_type,
                'material':             elem.material,
                'cylinder_length':      elem.cylinder_length,
                'cylinder_radius':      elem.cylinder_radius,
                'wedge_opening_angle':  elem.wedge_opening_angle,
                'wedge_height':         elem.wedge_height,
                'wedge_rotation_angle': elem.wedge_rotation_angle,
                'wedge_offset_x':       elem.wedge_offset_x,
                'wedge_offset_y':       elem.wedge_offset_y,
                'wedge_apex_to_base':   elem.wedge_apex_to_base,
                'fixed':                elem.fixed,
            })
        return pd.DataFrame(rows)

    def build_dipole_dataframe(self, global_z_offset: float = 0.0) -> pd.DataFrame:
        """Return a DataFrame with one row per Dipole."""
        rows = []
        for elem in self.elements:
            if not isinstance(elem, Dipole):
                continue
            rows.append({
                'name':             elem.name,
                'z_local':          elem.z_center,
                'z_global':         elem.z_center + global_z_offset,
                'field_strength':   elem.field_strength,
                'aperture':         elem.aperture,
                'length_z':         elem.length_z,
                'enge_coefficient': elem.enge_coefficient,
                'fixed':            elem.fixed,
            })
        return pd.DataFrame(rows)

    def build_elements_dataframe(self, global_z_offset: float = 0.0) -> pd.DataFrame:
        """Return a DataFrame with one row per element of any type."""
        rows = []
        for elem in self.elements:
            rows.append({
                'name':     elem.name,
                'type':     type(elem).__name__,
                'z_local':  elem.z_center,
                'z_global': elem.z_center + global_z_offset,
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
        polarities: Optional[List[int]] = None,
    ):
        super().__init__(name)
        self.n_cells = n_cells
        self.cell_length = cell_length
        self.total_length = total_length
        self.total_width = total_width
        self.polarities = polarities if polarities is not None else [1, -1, -1, 1]

        # Rendering / BDSIM field settings
        self.magnetic_field_model = magnetic_field_model
        self.magnetic_field_method = magnetic_field_method
        self.dipole_field_model = dipole_field_model
        self.interpolator = interpolator
        self.electric_field_model = electric_field_model
        self.reference_momentum = reference_momentum
        self.on_axis_tolerance = on_axis_tolerance

    def compute_rf_time_offsets(self, momentum_MeV_c: float, beam_start: float,
                                mass_MeV_c: float = 105.66) -> Dict[str, float]:
        """Compute and set time_offset on every RFCavity from transit time.

        t_i = (z_i_global - beam_start) / (beta * c)

        Parameters
        ----------
        momentum_MeV_c : float - Reference particle momentum [MeV/c]
        beam_start     : float - Global z position where the beam enters [m]
        mass_MeV_c     : float - Particle rest mass [MeV/c²] (default: muon)

        Returns
        -------
        dict mapping cavity name -> time_offset [ns]
        """
        C_M_PER_NS = 0.299792458  # speed of light [m/ns]

        E    = math.sqrt(momentum_MeV_c**2 + mass_MeV_c**2)
        beta = momentum_MeV_c / E

        global_offset = self.total_length / 2

        offsets = {}
        for elem in self.elements:
            if isinstance(elem, RFCavity):
                elem.time_offset = (elem.get_z(global_offset) - beam_start) / (beta * C_M_PER_NS)
                offsets[elem.name] = elem.time_offset

        return offsets

    def getCellStarts(self) -> List[List]:
        """Return the global z position and polarity sign at the start of each cell.

        Cell starts are evenly spaced across the full channel length, matching
        absorber placement positions, offset to global coordinates via total_length/2.
        The polarity sign ('+' or '-') is derived from self.polarities using the same
        stride logic as build_absorber_beamline, independent of whether absorbers or
        dipoles are actually present.

        Returns
        -------
        list of [z_global, sign] where sign is '+' or '-'
        """
        z_start = -self.n_cells * self.cell_length / 2
        global_offset = self.total_length / 2
        stride = max(1, len(self.polarities) // self.n_cells)
        result = []
        for cell_idx in range(self.n_cells):
            z = round(z_start + cell_idx * self.cell_length + global_offset, 2)
            pol = self.polarities[(cell_idx * stride) % len(self.polarities)]
            result.append([z, '+' if pol > 0 else '-'])
        return result

    def __repr__(self):
        return (f"MuonCoolingChannel(name='{self.name}', n_cells={self.n_cells}, "
                f"cell_length={self.cell_length}m, total_length={self.total_length}m, "
                f"total_width={self.total_width}m, n_elements={len(self.elements)})")
