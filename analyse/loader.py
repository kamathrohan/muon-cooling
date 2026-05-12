import argparse
from pathlib import Path

import pandas as pd
import numpy as np


def find_files(directory, pattern="*.txt"):
    """Recursively find all files matching pattern under directory."""
    return list(Path(directory).rglob(pattern))


def load_sims(filepaths):
    """
    Read each file and tag rows with SimName from the filename stem.
    Returns one concatenated DataFrame.
    """
    frames = []
    for path in filepaths:
        df = pd.read_csv(path)
        df["SimName"] = Path(path).stem
        frames.append(df)
    if not frames:
        raise ValueError("No files loaded")
    return pd.concat(frames, ignore_index=True)


def analyze_root_file(root_path: str, save=True):
    import pybdsim
    d = pybdsim.Data.Load(root_path)
    data_list = []
    for i in range(1, 200, 2):
        print(f"Processing sampler s{i}...")
        samplerdata = pybdsim.Data.SamplerData(d, 's' + str(i) + '.').data
        unique_s = np.unique(samplerdata['S'][samplerdata['S'] != 0.0])
        df = pd.DataFrame({
            "x":       samplerdata['x'],
            "y":       samplerdata['y'],
            "s":       unique_s[0],
            "indx":    i,
            "trackID": samplerdata['trackID'],
            "t":       samplerdata['T'],
            "E":       samplerdata['energy'],
            "p":       samplerdata['p'],
            "xp":      samplerdata['xp'],
            "yp":      samplerdata['yp'],
            "zp":      samplerdata['zp'],
        })
        data_list.append(df)
    final_df = pd.concat(data_list, ignore_index=True)
    if save:
        output_path = root_path.replace('.root', '.csv')
        final_df.to_csv(output_path, index=False)
        return final_df, output_path
    return final_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read a BDSIM ROOT file and write sampler data to CSV')
    parser.add_argument('--root', required=True, help='Path to input ROOT file')
    args = parser.parse_args()

    analyze_root_file(args.root)
