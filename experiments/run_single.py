
import argparse
from src.utils.config_loader import load_config
from src.simulation.run_simulation import run_simulation
import json

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
    simulation_data=run_simulation(cfg)
    with open(rf"C:\Users\Albi\Desktop\DnaA_manuscript\results\simulation_data.json", "w") as f:
        json.dump(simulation_data, f, indent=2)

    print("Loaded config:", args.config)
    print("LICENSING =", cfg.model.LICENSING)

if __name__ == "__main__":
    main()
