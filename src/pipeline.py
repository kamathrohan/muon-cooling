import argparse
import json
from itertools import groupby

from . import (
    MuonCoolingChannel,
    build_coil_beamline,
    build_dipole_beamline,
    build_absorber_beamline,
    build_rf_beamline,
    render_gmad,
)

_COIL_TEMPLATE_KEYS = {"name", "r_in", "r_thick", "N_pancakes", "L_pancake",
                       "L_spacing", "currDensity", "N_sheets", "material"}


def build_channel_from_config(config: dict) -> None:
    ch_cfg = config["channel"]
    channel = MuonCoolingChannel(
        n_cells=ch_cfg["n_cells"],
        cell_length=ch_cfg["cell_length"],
        total_length=ch_cfg["total_length"],
        total_width=ch_cfg["total_width"],
        reference_momentum=ch_cfg["reference_momentum"],
        on_axis_tolerance=ch_cfg.get("on_axis_tolerance", 2e-5),
    )

    for fixed_val, group in groupby(config.get("coils") or [], key=lambda c: c["fixed"]):
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
        channel = build_absorber_beamline(
            absorber_template=absorber_cfg,
            channel=channel,
        )

    rf_cfg = config.get("rf")
    if rf_cfg:
        rf_template = {k: v for k, v in rf_cfg.items() if k not in {"n_rf_cells", "rf_spacing"}}
        channel = build_rf_beamline(
            n_rf_cells=rf_cfg["n_rf_cells"],
            rf_spacing=rf_cfg["rf_spacing"],
            rf_template=rf_template,
            channel=channel,
        )

    out_cfg = config["output"]
    samp_cfg = config["samplers"]
    beam_cfg = config["beam"]

    render_gmad(
        channel,
        tpl_path=out_cfg["template_path"],
        out_gmad=out_cfg["gmad_path"],
        n_samplers=samp_cfg["n_samplers"],
        sampler_start_m=samp_cfg["sampler_start_m"],
        sampler_end_m=samp_cfg["sampler_end_m"],
        beam_mode=beam_cfg["mode"],
        beam_kwargs={"distr_file": beam_cfg["distr_file"]},
    )


def main():
    parser = argparse.ArgumentParser(description="Build a BDSIM channel .gmad from a JSON config")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    build_channel_from_config(config)


if __name__ == "__main__":
    main()
