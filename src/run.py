import argparse
import json
import os

import pybdsim

from .pipeline import build_channel_from_config


def run_from_config(config: dict) -> None:
    build_channel_from_config(config)

    out_cfg = config["output"]
    gmad_path = out_cfg["gmad_path"]
    output_name = out_cfg["output_name"]
    n_events = out_cfg["n_events"]

    pybdsim.Run.Bdsim(gmad_path, output_name, n_events)

    root_path = output_name + ".root"
    output_csv = config["analysis"]["output_csv"]
    n_particles = config["samplers"]["n_particles_per_sampler"]
    n_samplers = config["samplers"]["n_samplers"]

    _analyse(root_path, output_csv, n_samplers, n_particles)


def _analyse(root_path: str, output_csv: str, n_samplers: int, n_particles_per_sampler: int) -> None:
    import pandas as pd

    d = pybdsim.Data.Load(root_path)
    rows = []
    for i in range(1, n_samplers + 1):
        samplerdata = pybdsim.Data.SamplerData(d, f"s{i}.").data
        for j in range(n_particles_per_sampler):
            rows.append({
                "x": samplerdata["x"][j],
                "y": samplerdata["y"][j],
                "S": samplerdata["S"][j],
                "t": samplerdata["trackID"][j],
                "sampler_id": i,
            })

    pd.DataFrame(rows).to_csv(output_csv, index=False)


def main():
    parser = argparse.ArgumentParser(description="Build, run, and analyse a BDSIM channel from a JSON config")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    run_from_config(config)


if __name__ == "__main__":
    main()
