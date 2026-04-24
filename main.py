import argparse
from Models.Direct_SEIR import SEIRModel
from Models.Lattice_SEIR import Lattice ##importing both models from the Models folder


def main():
    parser = argparse.ArgumentParser(description="SEIR Simulation")   ##creates an arguement passer object allowing control over parameters
    parser.add_argument("--sigma", type=float, default=0.1)
    parser.add_argument("--gamma", type=float, default=0.005)
    parser.add_argument("--size", type=int, default=100)
    parser.add_argument("--total_particles", type=int, default=250)
    parser.add_argument("--ratio_exp", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=2000) ##control parameters added with default values included

    args = parser.parse_args() ##parse the arguements provided 

    # Deterministic model
    model = SEIRModel() ##create instance of the deterministic SEIR model
    model.calculate()
    model.plot() ##runs the functions

    # Lattice model
    lattice = Lattice( 
        size=args.size,
        sigma=args.sigma,
        gamma=args.gamma
    ) ##creates lattice object using the parameters from argparse

    lattice.initialise_grid(
        total_particles=args.total_particles,
        ratio_exp=args.ratio_exp
    ) 

    lattice.run(args.steps)
    lattice.plot_population()


if __name__ == "__main__":
    main() ##ensure main() only runs when this file is executed