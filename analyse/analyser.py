import argparse
import pybdsim
import pandas as pd


def analyze_root_file(root_path: str, output_csv_path: str):
    d = pybdsim.Data.Load(root_path)

    data_list = []
    for i in range(1, 120):
        samplerdata = pybdsim.Data.SamplerData(d, 's' + str(i) + '.').data
        for j in range(3):
            data_list.append({
                'x': samplerdata["x"][j],
                'y': samplerdata["y"][j],
                'S': samplerdata["S"][j],
                "t": samplerdata["trackID"][j],
                'sampler_id': i
            })

    df = pd.DataFrame(data_list)
    df.to_csv(output_csv_path, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze ROOT file and output sampler data to CSV')
    parser.add_argument('--root', required=True, help='Path to input ROOT file')
    parser.add_argument('--output', required=True, help='Path to output CSV file')
    args = parser.parse_args()

    analyze_root_file(args.root, args.output)
