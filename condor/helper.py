import os
import shutil
from typing import List

from jinja2 import Environment, FileSystemLoader

def makeFolders(simPath, jobID, iterNum, batchNum, nReplicas):
    """
    Create the necessary folders for storing simulation results.

    Parameters:
        simPath (str): Base path for simulations.
        jobID (str): Job ID.
        iterNum (str or int): Iteration number.
        batchNum (str or int): Batch (point) number.
        nReplicas (int): Number of replicas per batch point.

    Returns:
        List[str]: Paths to each replica folder.
    """
    simFolder = os.path.join(simPath, jobID)
    os.makedirs(simFolder, exist_ok=True)

    iterFolder = os.path.join(simFolder, str(iterNum))
    os.makedirs(iterFolder, exist_ok=True)

    if batchNum == 0:
        parent = iterFolder
    else:
        parent = os.path.join(iterFolder, str(batchNum))
        os.makedirs(parent, exist_ok=True)

    replicaFolders = []
    for r in range(nReplicas):
        replicaFolder = os.path.join(parent, str(r)+'r')
        os.makedirs(replicaFolder, exist_ok=True)
        replicaFolders.append(replicaFolder)
        jobOutputFolder = os.path.join(replicaFolder, "jobOutput")
        os.makedirs(jobOutputFolder, exist_ok=True)

    return replicaFolders


def renderDataGeneration(
    env_setup: str,
    bdsim_setup: str,
    out_path: str,
    tpl_path: str = "config/dataGeneration.sh.tpl",
) -> None:
    """
    Render the dataGeneration.sh script from its template.

    Parameters:
        env_setup (str): Path to the environment setup script (sourced first).
        bdsim_setup (str): Path to the BDSIM setup script (sourced second).
        out_path (str): Destination path for the rendered script.
        tpl_path (str): Path to the Jinja2 template file.
    """
    tpl_dir = os.path.dirname(os.path.abspath(tpl_path))
    tpl_file = os.path.basename(tpl_path)
    env = Environment(loader=FileSystemLoader(tpl_dir), keep_trailing_newline=True)
    script = env.get_template(tpl_file).render(
        env_setup=env_setup, bdsim_setup=bdsim_setup
    )
    with open(out_path, "w") as f:
        f.write(script)
    os.chmod(out_path, 0o755)


def renderSubmitArgs(
    max_runtime: int,
    out_path: str,
    tpl_path: str = "config/submitArgs.job.tpl",
) -> None:
    """
    Render the Condor submit file from its template.

    Parameters:
        max_runtime (int): Maximum job runtime in seconds (+MaxRuntime).
        out_path (str): Destination path for the rendered submit file.
        tpl_path (str): Path to the Jinja2 template file.
    """
    tpl_dir = os.path.dirname(os.path.abspath(tpl_path))
    tpl_file = os.path.basename(tpl_path)
    env = Environment(loader=FileSystemLoader(tpl_dir), keep_trailing_newline=True)
    script = env.get_template(tpl_file).render(max_runtime=max_runtime)
    with open(out_path, "w") as f:
        f.write(script)


def renderDoDataGeneration(
    replica_folders: List[str],
    jobcard: str,
    infile: str,
    ngenerate: int,
    out_path: str,
    tpl_path: str = "config/doDataGeneration.sh.tpl",
) -> None:
    """
    Render the doDataGeneration.sh script with hardcoded replica folder paths.

    Parameters:
        replica_folders (List[str]): Paths to each replica simulation folder.
        jobcard (str): Path to the rendered submitArgs.job file.
        infile (str): Input file name passed to bdsim via environment.
        ngenerate (int): Number of events to generate, passed via environment.
        out_path (str): Destination path for the rendered script.
        tpl_path (str): Path to the Jinja2 template file.
    """
    tpl_dir = os.path.dirname(os.path.abspath(tpl_path))
    tpl_file = os.path.basename(tpl_path)
    env = Environment(loader=FileSystemLoader(tpl_dir), keep_trailing_newline=True)
    script = env.get_template(tpl_file).render(
        replica_folders=replica_folders,
        jobcard=jobcard,
        infile=infile,
        ngenerate=ngenerate,
    )
    with open(out_path, "w") as f:
        f.write(script)
    os.chmod(out_path, 0o755)


def copyToFolders(src: str, folders: List[str]) -> None:
    """
    Copy a file into each folder in the list.

    Parameters:
        src (str): Path to the source file.
        folders (List[str]): Destination folders to copy the file into.
    """
    for folder in folders:
        shutil.copy2(src, folder)