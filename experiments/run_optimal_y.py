
import argparse
from src.utils.config_loader import load_config
from src.simulation.run_simulation import run_simulation
from src.utils.setup import get_n_star_n_forks, get_v_opt, get_chi0, get_alpha_opt, get_y_opt
import json
from dataclasses import replace


def make_optimal(cfg):
    n_star, _ = get_n_star_n_forks(cfg) 
    v_opt = get_v_opt(cfg, n_star)
    chi_0 = get_chi0(cfg, v_opt)
    cfg = replace(cfg, model=replace(cfg.model, CHI0=chi_0))
    alpha_opt=get_alpha_opt(cfg, v_opt)
    y_opt=get_y_opt(cfg, alpha_opt)
    cfg = replace(cfg, model=replace(cfg.model, COOP=y_opt))
    return cfg

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
    args = parser.parse_args()

    cfg = load_config(args.config)
    cfg = make_optimal(cfg)
    simulation_data=run_simulation(cfg)
    with open(rf"C:\Users\Albi\Desktop\DnaA_manuscript\results\simulation_data_optimal.json", "w") as f:
        json.dump(simulation_data, f, indent=2)

    print("Loaded config:", args.config)
    print("LICENSING =", cfg.model.LICENSING)

if __name__ == "__main__":
    main()
