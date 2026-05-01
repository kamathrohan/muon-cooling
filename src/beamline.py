import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Optional

from .elements import BeamlineElement, SolenoidCoil, RFCavity


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
