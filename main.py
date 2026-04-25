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
    parser.add_argument("--ratio_rec", type=float, default=0.0) 
    parser.add_argument("--steps", type=int, default=2000) ##control parameters added with default values included
    parser.add_argument("--lockdown_start", type=int, default= (10**100))
    parser.add_argument("--lockdown_strength", type=float, default = 0.0)
    parser.add_argument("--plot_step", type=int, default=0)

    parser.add_argument("--exposure_rate", type=float, default = 1.0)
    parser.add_argument("--incubation_rate", type=float, default = 1.0)
    parser.add_argument("--recovery_rate", type=float, default = 0.1)

    parser.add_argument("--sus0", type=float, default=0.99)
    parser.add_argument("--exp0", type=float, default=0.01)
    parser.add_argument("--inf0", type=float, default=0.0)
    parser.add_argument("--rec0", type=float, default=0.0)

    parser.add_argument("--days", type=int, default=100) ##parameters of the determinsitic model added


    args = parser.parse_args() ##parse the arguements provided 

    # Deterministic model
    model = SEIRModel(
        exposure_rate=args.exposure_rate,
        incubation_rate=args.incubation_rate,
        recovery_rate=args.recovery_rate,
        sus0=args.sus0,
        exp0=args.exp0,
        inf0=args.inf0,
        rec0=args.rec0,
        days=args.days
    ) ##create instance of the deterministic SEIR model
    model.calculate()
    model.plot() ##runs the functions

    # Lattice model
    lattice = Lattice( 
        size=args.size,
        sigma=args.sigma,
        gamma=args.gamma,
        lockdown_start=args.lockdown_start,
        lockdown_strength=args.lockdown_strength

    ) ##creates lattice object using the parameters from argparse

    lattice.initialise_grid(
        total_particles=args.total_particles,
        ratio_exp=args.ratio_exp,
        ratio_rec=args.ratio_rec
    ) 

    lattice.run(args.steps, args.plot_step)
    lattice.plot_population()



if __name__ == "__main__":
    main() ##ensure main() only runs when this file is executed