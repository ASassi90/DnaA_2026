
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


def get_range(sweep_dict):
    if sweep_dict["scale"]=="log":
        maxexp=np.log(sweep_dict["max_val"])
        minexp=np.log(sweep_dict["min_val"])
        steps=sweep_dict["steps"]
        values=np.exp(np.linspace(minexp, maxexp, steps))
    else:
        max_v=sweep_dict["max_val"]
        min_v=sweep_dict["min_val"]
        steps=sweep_dict["steps"]
        values=np.linspace(min_v, max_v, steps)    
    
    return values

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
        help="Path to the YAML sweep"
    )
    args = parser.parse_args()
    sweep_dict=yaml.safe_load(open(args.sweep))
    print("base path = ",sweep_dict["base_yaml"])
    params=sweep_dict["params"]
    sweep_coop=params["COOP"]
    sweep_change=params["CHANGE"]
    y_values=get_range(sweep_coop)
    change_values=get_range(sweep_change)
    for change in change_values:
        cfg = load_config(args.config)
        cfg = replace(cfg, model=replace(cfg.model, CHANGE=change))
        cfg = make_optimal(cfg)
        for y_new in y_values:
            print("cooperativity strength (y): ", y_new)
            kori_new=change_kori(cfg, y_new)
            cfg0=replace(cfg, model=replace(cfg.model, COOP=y_new, K_OPEN=kori_new))
            simulation_data=run_simulation(cfg0)
            with open(rf"C:\Users\Albi\Desktop\DnaA_manuscript\results\data\optimal_y_%g_change_%g.json"%(y_new, cfg.model.CHANGE), "w") as f:
                json.dump(simulation_data, f, indent=2)
            print("chi0 and coop = ", cfg0.model.CHI0, cfg0.model.COOP)
    print("Loaded config:", args.config)
    print("LICENSING =", cfg.model.LICENSING)

if __name__ == "__main__":
    main()
