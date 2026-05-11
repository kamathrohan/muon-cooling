#!/usr/bin/env python3
import json
import os
import shutil
import sys
from typing import List

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.pipeline import build_channel_from_config
from condor.helper import makeFolders, renderDataGeneration, renderSubmitArgs, renderDoDataGeneration


def setup_simulation(simDir: str, channel_json: str, template: str, beam_file: str,
                     gmad_name: str = None, distr_file_override: str = None,
                     extra_config: dict = None) -> str:
    """
    Create simDir, render a .gmad file into it, and copy the beam distribution file.

    Parameters:
        simDir (str): Directory to create and write the simulation into.
        channel_json (str): Path to the channel JSON config file.
        template (str): Path to the Jinja2 .tpl template file.
        beam_file (str): Path to the beam distribution file to copy into simDir.
        gmad_name (str): Override for the output .gmad filename (without directory).
        distr_file_override (str): If provided, use this path in the gmad distrFile instead
            of the one in the config, and skip copying the beam file into simDir.
        extra_config (dict): Key/value pairs merged into the loaded config before rendering.
            Nested dicts are merged shallowly (one level deep).

    Returns:
        str: Path to the rendered .gmad file.
    """
    os.makedirs(simDir, exist_ok=True)

    with open(channel_json) as f:
        config = json.load(f)

    if extra_config:
        for k, v in extra_config.items():
            if isinstance(v, dict) and isinstance(config.get(k), dict):
                config[k].update(v)
            else:
                config[k] = v

    expected_beam = config["beam"]["distr_file"]
    if os.path.basename(beam_file) != os.path.basename(expected_beam):
        raise ValueError(
            f"beam_file '{os.path.basename(beam_file)}' does not match "
            f"distr_file '{expected_beam}' in {channel_json}"
        )

    if gmad_name is None:
        gmad_name = os.path.splitext(os.path.basename(channel_json))[0] + ".gmad"

    gmad_path = os.path.join(simDir, gmad_name)

    if distr_file_override is not None:
        config["beam"]["distr_file"] = distr_file_override
    else:
        shutil.copy2(beam_file, simDir)

    build_channel_from_config(config, template, gmad_path)

    return gmad_path


def initialiseJob(
    simDir: str,
    outputDir: str,
    jobId: str,
    iterNum,
    nBatch: int,
    nReplicas: int,
    channel_json: str,
    beam_file: str,
    channel_template: str,
    datagen_template: str,
    submit_template: str,
    do_datagen_template: str,
    max_runtime: int,
    env_setup: str,
    bdsim_setup: str,
    ngenerate: int,
) -> List[str]:
    """
    Set up all batch/replica simulation folders and render the dataGeneration script.

    Parameters:
        simDir (str): Base simulation directory; job folder created inside it.
        outputDir (str): Base output directory; job folder created inside it.
        jobId (str): Job identifier used as the top-level folder name.
        iterNum: Iteration number.
        nBatch (int): Number of independent perturbation draws (batches).
        nReplicas (int): Number of replica folders per batch (share the same draw).
        channel_json (str): Path to the channel JSON config file.
        beam_file (str): Path to the beam distribution file.
        channel_template (str): Path to the Jinja2 channel .tpl template.
        datagen_template (str): Path to the dataGeneration.sh.tpl template.
        submit_template (str): Path to the submitArgs.job.tpl template.
        do_datagen_template (str): Path to the doDataGeneration.sh.tpl template.
        max_runtime (int): Maximum job runtime in seconds (+MaxRuntime).
        env_setup (str): Path to the environment setup script.
        bdsim_setup (str): Path to the BDSIM setup script.
        ngenerate (int): Number of events exported for bdsim.

    Returns:
        List[str]: Paths to all replica folders across all batches.
    """
    beam_basename = os.path.basename(beam_file)
    beam_rel = "../" + beam_basename

    rng = np.random.default_rng()
    all_folders = []
    all_stems = []

    for b in range(nBatch):
        batch_folders = makeFolders(simDir, jobId, iterNum, b, nReplicas)
        batch_stems = [
            f"channel_{jobId}_{iterNum}n_{b}b_{r}r"
            for r in range(nReplicas)
        ]

        beam_parent = os.path.dirname(batch_folders[0])
        shutil.copy2(beam_file, beam_parent)

        tol_seed = int(rng.integers(0, 2**31))
        extra = {"toleranceCoil": {"seed": tol_seed}}

        for folder, stem in zip(batch_folders, batch_stems):
            setup_simulation(folder, channel_json, channel_template, beam_file,
                             gmad_name=stem + ".gmad", distr_file_override=beam_rel,
                             extra_config=extra)

        all_folders.extend(batch_folders)
        all_stems.extend(batch_stems)

    job_output_dir = os.path.join(outputDir, jobId)
    os.makedirs(job_output_dir, exist_ok=True)

    datagen_out = os.path.join(job_output_dir, "dataGeneration.sh")
    renderDataGeneration(
        env_setup=env_setup,
        bdsim_setup=bdsim_setup,
        out_path=datagen_out,
        tpl_path=datagen_template,
    )

    jobcard = os.path.join(job_output_dir, "submitArgs.job")
    renderSubmitArgs(
        max_runtime=max_runtime,
        out_path=jobcard,
        tpl_path=submit_template,
    )

    renderDoDataGeneration(
        replica_folders=all_folders,
        outfiles=all_stems,
        jobcard=jobcard,
        ngenerate=ngenerate,
        out_path=os.path.join(job_output_dir, "doDataGeneration.sh"),
        tpl_path=do_datagen_template,
    )

    return all_folders


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Initialise all replica folders and render the dataGeneration script")
    parser.add_argument("--config", required=True, help="Path to job JSON config file")
    parser.add_argument("--job-id", required=True, help="Job identifier")
    parser.add_argument("--run", action="store_true", default=False, help="Execute doDataGeneration.sh after setup")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    replica_folders = initialiseJob(
        simDir=cfg["simDir"],
        outputDir=cfg["outputDir"],
        jobId=args.job_id,
        iterNum=cfg["iterNum"],
        nBatch=cfg["nBatch"],
        nReplicas=cfg["nReplicas"],
        channel_json=cfg["channel_json"],
        beam_file=cfg["beam_file"],
        channel_template=cfg["channel_template"],
        datagen_template=cfg["datagen_template"],
        submit_template=cfg["submit_template"],
        do_datagen_template=cfg["do_datagen_template"],
        max_runtime=cfg["max_runtime"],
        env_setup=cfg["env_setup"],
        bdsim_setup=cfg["bdsim_setup"],
        ngenerate=cfg["ngenerate"],
    )
    print(f"Job '{args.job_id}' ready with {len(replica_folders)} replicas:")
    for folder in replica_folders:
        print(f"  {folder}")

    if args.run:
        job_output_dir = os.path.join(os.path.abspath(cfg["outputDir"]), args.job_id)
        do_datagen = os.path.join(job_output_dir, "doDataGeneration.sh")
        os.chdir(job_output_dir)
        os.execv("/bin/bash", ["/bin/bash", do_datagen])


if __name__ == "__main__":
    main()
