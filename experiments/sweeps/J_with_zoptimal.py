
import argparse
from src.utils.config_loader import load_config
from src.simulation.run_simulation import run_simulation
from experiments.run_optimal_y import make_optimal
import json
from dataclasses import replace
import yaml
import numpy as np

def change_kori(cfg, y):
    """
        This function returns a new value of K_ori, still optimal, corresponding to a new 
        value of the cooperativity strength y, all other parameters being kept the same. 
    """
    epsilon_cost, origin_sites, K, dnaa, kori_old, y_old =  (cfg.model.E_COST,
                                                             cfg.model.ORIGIN_SITES,
                                                             cfg.model.K,
                                                             cfg.model.DNAA_CONCENTRATION,
                                                             cfg.model.K_OPEN,
                                                             cfg.model.COOP)
    alpha=kori_old*np.exp(epsilon_cost/origin_sites)/(y_old*np.sqrt(K*dnaa))
    kori =y/(np.exp(epsilon_cost/origin_sites)/(alpha*np.sqrt(K*dnaa)))
    return kori


def main():
    parser = argparse.ArgumentParser(
        description="Run one simulation from a YAML configuration."
    )
    parser.add_argument(
        "--config",
        required=True,
        metavar="YAML_FILE",
        help="Path to the YAML config (e.g., configs/base.yaml)"
    )
    parser.add_argument(
        "--sweep",
        required=True,
        metavar="YAML_FILE",
        help="Path to the YAML sweep (e.g., configs/base.yaml)"
    )
    args = parser.parse_args()
    sweep_dict=yaml.safe_load(open(args.sweep))
    if sweep_dict["scale"]=="log":
        maxexp=np.log(sweep_dict["max_val"])
        minexp=np.log(sweep_dict["min_val"])
        steps=sweep_dict["steps"]
        y_values=np.exp(np.linspace(minexp, maxexp, steps))
    else:
        max_v=sweep_dict["max_val"]
        min_v=sweep_dict["min_val"]
        steps=sweep_dict["steps"]
        y_values=np.exp(np.linspace(min_v, max_v, steps))        

    cfg = load_config(args.config)
    cfg = replace(cfg, model=replace(cfg.model, CHANGE=1.05))
    cfg = make_optimal(cfg)
    for y_new in y_values:
        print("cooperativity strength (y): ", y_new)
        kori_new=change_kori(cfg, y_new)
        cfg0=replace(cfg, model=replace(cfg.model, COOP=y_new, K_OPEN=kori_new))
        simulation_data=run_simulation(cfg0)
        with open(rf"C:\Users\Albi\Desktop\DnaA_manuscript\results\data\optimal_y_%g_chi0_%g.json"%(y_new, cfg.model.CHI0), "w") as f:
            json.dump(simulation_data, f, indent=2)

    print("Loaded config:", args.config)
    print("LICENSING =", cfg.model.LICENSING)

if __name__ == "__main__":
    main()
