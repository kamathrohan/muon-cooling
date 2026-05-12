import argparse
import json
from itertools import groupby

import numpy as np

from . import (
    MuonCoolingChannel,
    build_coil_beamline,
    build_dipole_beamline,
    build_absorber_beamline,
    build_rf_beamline,
    render_gmad,
)
from .physics.elements import SolenoidCoil

_COIL_TEMPLATE_KEYS = {"name", "r_in", "r_thick", "N_pancakes", "L_pancake",
                       "L_spacing", "currDensity", "N_sheets", "material"}

_ABSORBER_BUILDER_KEYS = {"n_cells", "wedge_alignment_angle1", "wedge_alignment_angle2",
                          "wedge_offset_x"}

_RF_BUILDER_KEYS = {"n_cells", "n_rf_cells", "rf_spacing"}


def build_channel_from_config(config: dict, tpl_path: str, out_gmad: str) -> None:
    ch_cfg = config["channel"]
    channel = MuonCoolingChannel(
        n_cells=ch_cfg["n_cells"],
        cell_length=ch_cfg["cell_length"],
        total_length=ch_cfg["total_length"],
        total_width=ch_cfg["total_width"],
        reference_momentum=ch_cfg["reference_momentum"],
        on_axis_tolerance=ch_cfg.get("on_axis_tolerance", 2e-5),
        magnetic_field_model=ch_cfg.get("magnetic_field_model", "solenoidsheet"),
        magnetic_field_method=ch_cfg.get("magnetic_field_method", "grid"),
        grid_points_per_mm=ch_cfg.get("grid_points_per_mm", 0.1),
        z_period_start=ch_cfg.get("z_period_start", -70000),
        z_period_end=ch_cfg.get("z_period_end", 175000),
        period_length=ch_cfg.get("period_length", 2000),
    )

    for fixed_val, group in groupby(sorted(config.get("coils") or [], key=lambda c: c["fixed"]), key=lambda c: c["fixed"]):
        group_list = list(group)
        channel = build_coil_beamline(
            z_coils=[c["z_position"] for c in group_list],
            coil_templates=[{k: v for k, v in c.items() if k in _COIL_TEMPLATE_KEYS}
                            for c in group_list],
            polarities=[c["polarity"] for c in group_list],
            fixed=fixed_val,
            channel=channel,
        )

    dipole_cfg = config.get("dipole")
    if dipole_cfg:
        coil_cell_z = dipole_cfg["coil_cell_z"]
        dipole_template = {k: v for k, v in dipole_cfg.items() if k != "coil_cell_z"}
        channel = build_dipole_beamline(
            coil_cell_z=coil_cell_z,
            dipole_template=dipole_template,
            channel=channel,
        )

    absorber_cfg = config.get("absorber")
    if absorber_cfg:
        absorber_template = {k: v for k, v in absorber_cfg.items() if k not in _ABSORBER_BUILDER_KEYS}
        channel = build_absorber_beamline(
            n_cells=absorber_cfg.get("n_cells"),
            wedge_alignment_angle1=absorber_cfg.get("wedge_alignment_angle1", 3.141592653589793),
            wedge_alignment_angle2=absorber_cfg.get("wedge_alignment_angle2", 0.0),
            offsetX=absorber_cfg.get("wedge_offset_x", 0.0),
            absorber_template=absorber_template,
            channel=channel,
        )

    rf_cfg = config.get("rf")
    if rf_cfg:
        rf_template = {k: v for k, v in rf_cfg.items() if k not in _RF_BUILDER_KEYS}
        channel = build_rf_beamline(
            n_cells=rf_cfg.get("n_cells"),
            n_rf_cells=rf_cfg["n_rf_cells"],
            rf_spacing=rf_cfg["rf_spacing"],
            rf_template=rf_template,
            channel=channel,
        )

    rf_phasing = config.get("rf_phasing")
    rf_phases  = config.get("rf_phases")
    if rf_phasing and rf_phases:
        raise ValueError("Specify either 'rf_phasing' (computed) or 'rf_phases' (explicit list), not both.")
    if rf_phasing:
        mode = rf_phasing.get("mode")
        if mode == "compute":
            channel.compute_rf_time_offsets(
                momentum_MeV_c=rf_phasing["momentum_MeV_c"],
                beam_start=rf_phasing["beam_start"],
                mass_MeV_c=rf_phasing.get("mass_MeV_c", 105.66),
            )
        else:
            raise ValueError(f"Unknown rf_phasing mode '{mode}'. Only 'compute' is supported.")
    elif rf_phases:
        channel.set_rf_time_offsets(rf_phases)

    coil_tilts = config.get("coilTilts")
    if coil_tilts:
        n_coils = sum(1 for e in channel.elements if isinstance(e, SolenoidCoil))
        expanded = {}
        for key in ("tiltX", "tiltY", "tiltZ"):
            val = coil_tilts[key]
            if isinstance(val, list):
                if len(val) != n_coils:
                    raise ValueError(
                        f"coilTilts.{key} has {len(val)} values but channel has {n_coils} coils."
                    )
                expanded[key] = val
            else:
                expanded[key] = [val] * n_coils
        channel.set_tilts(expanded)

    tol_cfg = config.get("tolerance", {}).get("coil")
    if tol_cfg:
        rng = np.random.default_rng(tol_cfg.get("seed"))

        def _draw(sigma):
            return np.clip(rng.normal(0, sigma), -3 * sigma, 3 * sigma)

        for coil in (e for e in channel.elements if isinstance(e, SolenoidCoil)):
            if tol_cfg.get("current", 0):
                r = _draw(tol_cfg["current"])
                coil.currDensity *= (1 + r)
                coil._Itot_per_pancake *= (1 + r)
            if tol_cfg.get("tilt", 0):
                coil.tilt_x += _draw(tol_cfg["tilt"])
                coil.tilt_y += _draw(tol_cfg["tilt"])
            if tol_cfg.get("offset", 0):
                coil.z_center *= (1 + _draw(tol_cfg["offset"]))

    samp_cfg = config["samplers"]
    beam_cfg = config["beam"]

    render_gmad(
        channel,
        tpl_path=tpl_path,
        out_gmad=out_gmad,
        sampler_mode=samp_cfg.get("sampler_mode", "linspace"),
        sampler_kwargs=samp_cfg.get("sampler_kwargs"),
        beam_mode=beam_cfg["mode"],
        beam_kwargs={"distr_file": beam_cfg["distr_file"]},
    )


def main():
    parser = argparse.ArgumentParser(description="Build a BDSIM channel .gmad from a JSON config")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--template", required=True, help="Path to Jinja2 .tpl template file")
    parser.add_argument("--output", required=True, help="Output path for the rendered .gmad file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    build_channel_from_config(config, args.template, args.output)


if __name__ == "__main__":
    main()
