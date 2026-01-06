import yaml
import numpy as np
from experiments.sweeps.run_y_and_chi0 import get_range
import argparse
import json
from matplotlib import pyplot as plt
from typing import Sequence, Tuple, List
from src.utils.helpers import create_figure

def get_discontinuities(origins: Sequence[float], forks: Sequence[float]):                        
    """
    Detect indices of initiations, terminations, and divisions based on two traces:
    origins and forks. Indices returned correspond to the right-hand sample (i+1),
    i.e., the index where the new value appears after a change between i -> i+1.

    Rules:
      - Initiation:  Δorigins > 0 and Δforks > 0
      - Division:    Δorigins < 0 and Δforks < 0
      - Termination: Δorigins == 0 and Δforks < 0

    Parameters
    ----------
    origins : Sequence[float]
        Time series of origin counts.
    forks : Sequence[float]
        Time series of fork counts.

    Returns
    -------
    initiations, terminations, divisions : Tuple[List[int], List[int], List[int]]
        Lists of indices (0-based) where each event is detected. The index refers
        to the right-hand point of the transition (i+1).

    Raises
    ------
    ValueError
        If lengths differ or sequences are shorter than 2.

    Notes
    -----
    - If you prefer 1-based indices, add +1 to each reported index.
    - If the data is noisy or floating-point, consider adding a tolerance.
    """

    if len(origins) != len(forks):
        raise ValueError("origins and forks must have the same length.")
    if len(origins) < 2:
        raise ValueError("origins and forks must have length >= 2.")

    initiations: List[int] = []
    terminations: List[int] = []
    divisions: List[int] = []

    # Scan adjacent pairs (i -> i+1), record i+1 on matching rule
    for i in range(len(origins) - 1):
        d_o = origins[i + 1] - origins[i]
        d_f = forks[i + 1] - forks[i]

        if d_o > 0 and d_f > 0:
            initiations.append(i + 1)
        elif d_o < 0 and d_f < 0:
            divisions.append(i + 1)
        elif d_o == 0 and d_f < 0:
            terminations.append(i + 1)

    return initiations, terminations, divisions



def main():
    parser = argparse.ArgumentParser(
        description="Run one simulation from a YAML configuration."
    )
    parser.add_argument(
        "--sweep",
        required=True,
        metavar="YAML_FILE",
        help="Path to the YAML sweep"
    )
    args = parser.parse_args()
    sweep_dict=yaml.safe_load(open(args.sweep))
    params=sweep_dict["params"]
    sweep_coop=params["COOP"]
    sweep_change=params["CHANGE"]
    y_values=get_range(sweep_coop)
    change_values=get_range(sweep_change)
    plt.figure()
    for change in change_values:
        print(change)
        cv_volumes=[]
        for y_new in y_values:
            file_path=rf"C:\Users\Albi\Desktop\DnaA_manuscript\results\data\optimal_y_%g_change_%g.json"%(y_new, change)
            with open(file_path, "r") as f:
                simulation_data=json.load(f)
                #time=simulation_data["time"]
                volume=simulation_data["volume"]
                origins=simulation_data["origins"]
                forks=simulation_data["n_forks"]
                initiations, _, _ = get_discontinuities(origins, forks)
                volume_at_in=[volume[ind] for ind in initiations]
                vol=volume_at_in[int(len(volume_at_in)/3):]
                cv=np.sqrt(np.var(vol))/np.mean(vol)
                cv_volumes.append(cv)
                """
                plt.figure()
                plt.plot(time, origins)
                initiation_times=[time[in_index] for in_index in initiations]
                for t in initiation_times:
                    plt.axvline(t)

                plt.show()
                plt.close()
                """
        print("CV of initiation volume= ", cv_volumes)
        plt.plot(np.log(y_values), np.log(cv_volumes), '-o', label=rf"$\chi/V^*$ =%.2g"%(1./change))
        #plt.xscale("log")
        #plt.yscale("log")
        plt.legend()

    plt.show()

if __name__=="__main__":
    main()