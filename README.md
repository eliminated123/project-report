# SEIR Epidemic Simulation

A Python implementation of SEIR (Susceptible → Exposed → Infected → Recovered) epidemic modelling using two different approaches: a deterministic based model and a stochastic Monte Carlo lattice model.

## Overview

The SEIR framework divides a population into four compartments:

| State | Symbol | Description |
|-------|--------|-------------|
| Susceptible | S | Can be infected |
| Exposed | E | Infected but not yet infectious |
| Infected | I | Infectious and can spread the disease |
| Recovered | R | Recovered and immune |

This project implements the SEIR framework using two different models:

- **Deterministic model** — differential equations operating on population fractions.
- **Monte Carlo lattice model** — individual agents placed on a 2D grid with spatial movement and local interactions which captures stochastic and spatial effects.

---

## Models

### Deterministic SEIR Model

Located in `Models/Direct_SEIR.py`.

The model evolves population fractions using discrete difference equations at each timestep (days):


- S(t+1) = S(t) - β · I(t) · S(t)
- E(t+1) = E(t) + β · I(t) · S(t) - σ · E(t)
- I(t+1) = I(t) + σ · E(t) - γ · I(t)
- R(t+1) = R(t) + γ · I(t)


Where:
- **β** (`infection_rate`) —  rate of transmission (rate at which susceptibles become exposed)
- **σ** (`incubation_rate`) — rate at which exposed individuals become infectious 
- **γ** (`recovery_rate`) — rate at which infected individuals recover

Population conservation (`S + E + I + R = 1`) is enforced at every step.

---

### Monte Carlo Lattice Model

Located in `Models/Lattice_SEIR.py`.

Agents are placed on a 2D `size × size` grid (hard wall boundary conditions). Each cell is either empty (0) or occupied by an agent in one of four states:

| Value | State |
|-------|-------|
| 0 | Empty |
| 1 | Susceptible |
| 2 | Exposed |
| 3 | Infected |
| 4 | Recovered |

Each Monte Carlo step consists of:

1. **Movement** — Every agent attempts to move to a randomly selected neighbouring cell (4 orthogonal + 4 diagonal) if it is vacant. Movement can be restricted by lockdown.
2. **State updates** — Agent states are updated based on local neighbourhood interactions:
   - **S → E**: Probability of updating state is proportional to the fraction of infected neighbours.
   - **E → I**: With probability `sigma` per step.
   - **I → R**: With probability `gamma` per step.
3. **Recording** — Population counts for each state are appended to their corresponding history lists.
4. **Validity check** — Population conservation and states are verified.

---

## Project Structure

```
project/
│
├── main.py                   # Entry point — runs both models with command line arguments
│
├── Models/
│   ├── Direct_SEIR.py        # Deterministic ODE-based SEIR model
│   └── Lattice_SEIR.py       # Monte Carlo agent-based lattice model
│
└── README.md
```

---

## Installation

**Requirements:** 

The following libraries:

```
numpy
matplotlib
argparse
```

Install:

```bash
pip install numpy matplotlib argparse
```

Clone the repository:

```bash
git clone 2588092.bundle project-report
cd project-report
```

---

## Usage

### Running via Command Line

Run both models with default parameters:

```bash
python main.py
```

Run with custom parameters:

```bash
python main.py --sigma 0.15 --gamma 0.01 --steps 1000 --total_particles 500
```

Simulate a lockdown starting at step 200 with 80% movement reduction:

```bash
python main.py --lockdown_start 200 --lockdown_strength 0.2
```

Simulate a vaccinated initial population (10% recovered at start):

```bash
python main.py --ratio_rec 0.1
```

Plot the lattice grid at a specific step:

```bash
python main.py --plot_step 100
```

---

### Using Models in Python

**Deterministic model:**

```python
from Models.Direct_SEIR import SEIRModel

model = SEIRModel(
        infection_rate=args.infection_rate,
        incubation_rate=args.incubation_rate,
        recovery_rate=args.recovery_rate,
        sus0=args.sus0,
        exp0=args.exp0,
        inf0=args.inf0,
        rec0=args.rec0,
        days=args.days
) 
model.calculate()
model.plot() 
```

**Lattice model:**

```python
lattice = Lattice(
    size=args.size,
    sigma=args.sigma,
    gamma=args.gamma,
    lockdown_start=args.lockdown_start,
    lockdown_strength=args.lockdown_strength
)
lattice.initialise_grid(
    total_particles=args.total_particles,
    ratio_exp=args.ratio_exp,
    ratio_rec=args.ratio_rec
)
lattice.run(args.steps, args.plot_step)
lattice.plot_population()
```

---

## Parameters

### Deterministic Model Parameters

| Parameter | CLI Flag | Default | Description |
|-----------|----------|---------|-------------|
| Infection rate (β) | `--infection_rate` | `1.0` | Transmission rate from infected → exposed  |
| Incubation rate (σ) | `--incubation_rate` | `1.0` | Rate of exposed → infected transition |
| Recovery rate (γ) | `--recovery_rate` | `0.1` | Rate of infected → recovered transition |
| Initial susceptible | `--sus0` | `0.99` | Starting fraction of susceptible population |
| Initial exposed | `--exp0` | `0.01` | Starting fraction of exposed population |
| Initial infected | `--inf0` | `0.0` | Starting fraction of infected population |
| Initial recovered | `--rec0` | `0.0` | Starting fraction of recovered population |
| Days | `--days` | `100` | Number of timesteps to simulate |

> **Note:** `sus0 + exp0 + inf0 + rec0` must equal `1.0`.

---

### Lattice Model Parameters

| Parameter | CLI Flag | Default | Description |
|-----------|----------|---------|-------------|
| Grid size | `--size` | `100` | Lattice dimensions (`size × size`) |
| Sigma (σ) | `--sigma` | `0.1` | Probability per step of E → I transition |
| Gamma (γ) | `--gamma` | `0.005` | Probability per step of I → R transition |
| Total particles | `--total_particles` | `250` | Number of agents placed on the grid |
| Exposed ratio | `--ratio_exp` | `0.05` | Fraction of agents initially exposed |
| Recovered ratio | `--ratio_rec` | `0.0` | Fraction of agents initially recovered (vaccinated) |
| MC steps | `--steps` | `2000` | Number of Monte Carlo steps |
| Plot step | `--plot_step` | `0` | Step at which to plot the grid (0 = no grid plot) |
| Lockdown start | `--lockdown_start` | `∞` | Step at which lockdown begins |
| Lockdown strength | `--lockdown_strength` | `0.0` | Movement probability during lockdown (0.0 = full lockdown, 1.0 = no effect) |

---

## Interventions

### Vaccination

Fraction of the population immunised by setting `ratio_rec > 0`. These agents begin in the recovered (immune) state and cannot be infected. Combined with `ratio_exp`, `ratio_exp + ratio_rec ≤ 1` must hold.

```bash
python main.py --ratio_rec 0.2   # 20% of agents vaccinated at start
```

### Lockdown

A time-dependent lockdown restricts agent movement after a chosen step. `lockdown_strength` defines the probability that each agent is permitted to move during the lockdown period:

- `lockdown_strength = 0.0` → full lockdown (no movement)
- `lockdown_strength = 1.0` → no restriction 
- `lockdown_strength = 0.5` → 50% movement probability

```bash
python main.py --lockdown_start 300 --lockdown_strength 0.5
```

---

## Output & Visualisation

### Deterministic Model

`model.plot()` produces a line plot of population fractions (S, E, I, R) over time in days.

### Lattice Model

`lattice.plot_population()` produces a line plot of agent counts (S, E, I, R) over Monte Carlo steps.

`lattice.plot_grid()` (when using `--plot_step`) produces a 2D scatter plot of agent positions on the lattice, colour coded by state:

| Colour | State |
|--------|-------|
| Blue | Susceptible |
| Orange | Exposed |
| Red | Infected |
| Green | Recovered |

---

## Class Design & Encapsulation

The two classes use different but well motivated levels of encapsulation that reflect 
their internal complexity and intended public interfaces.

---

### `SEIRModel` — Strong Private Encapsulation

`SEIRModel` uses Python's double underscore convention (`__`) to make its internal 
data and update methods strictly private:

| Member | Access | Purpose |
|--------|--------|---------|
| `__init__()` | Public | Configure model parameters |
| `calculate()` | Public | Run the simulation |
| `plot()` | Public | Visualise results |
| `check_validity()` | Public | Verify conservation at each step |
| `__sus_list`, `__exp_list`, `__inf_list`, `__rec_list` | Private (`__`) | Internal population history — fully protected |
| `__update_sus/exp/inf/rec()` | Private (`__`) | Difference equation steps — internal only |

This is well-motivated because the population fraction lists must always satisfy 
S + E + I + R = 1.0. Allowing external code to modify them would risk breaking 
this conservation. The update methods are likewise hidden as they must only 
ever be called in a specific sequence inside `calculate()`. The result is a minimal 
public interface — a user instantiates the class, calls `calculate()`, then `plot()`.

---

### `Lattice` — Protected Internal Methods with Public Data

`Lattice` uses the single-underscore convention (`_`) for its internal methods, marking 
them as protected rather than fully private. This is a deliberate middle ground, the 
methods are signalled as internal implementation details not intended for direct external 
use, while remaining accessible if needed for subclassing or testing.

| Member | Access | Purpose |
|--------|--------|---------|
| `grid`, `sus_pop`, `exp_pop`, `inf_pop`, `rec_pop` | Public | Grid state and population history — inspectable externally |
| `sigma`, `gamma`, `lockdown_start`, etc. | Public | Tuneable simulation parameters |
| `initialise_grid()` | Public | Place agents on the grid |
| `step()` / `run()` | Public | Advance the simulation |
| `plot_population()` / `plot_grid()` | Public | Visualisation |
| `_move_agents()` | Protected (`_`) | Internal movement logic |
| `_update_agents()` | Protected (`_`) | Internal state transition logic |
| `_record_counts()` | Protected (`_`) | Internal population recording |
| `_check_validity()` | Protected (`_`) | Internal integrity check |
| `_get_neighbours()` | Protected (`_`) | Internal neighbourhood query |

The five protected methods are sub steps that must only be called in the correct order 
through `step()`. Marking them as protected communicates this without completely 
restricting access. The data attributes remain fully public as they can be 
inspected externally without any risk of effecting the simulation state.

---

### Summary

| | `SEIRModel` | `Lattice` |
|---|---|---|
| Internal data | Private (`__`) | Public |
| Internal methods | Private (`__`) | Protected (`_`) |


## Validity Checks

Both models perform validity checks at each step to ensure realistic results:

**Deterministic model** (`check_validity`):
- Total population fractions sum to `1.0` (within tolerance `1e-6`).
- No population fraction is negative.

**Lattice model** (`check_validity`):
- Total agent count equals `total_particles`.
- All grid values are within the valid range `[0, 4]`.
- `ratio_exp` and `ratio_rec` are each within `[0, 1]`.
- `ratio_exp + ratio_rec ≤ 1`.

A `ValueError` is raised immediately if any check fails.

---

## Examples

**Basic epidemic with no intervention:**
```bash
python main.py --sigma 0.1 --gamma 0.005 --steps 2000
```

**Epidemic with early lockdown:**
```bash
python main.py --lockdown_start 100 --lockdown_strength 0.05 --steps 2000
```

**High vaccination coverage:**
```bash
python main.py --ratio_rec 0.4 --ratio_exp 0.05 --steps 2000
```

**Snapshot of the spatial grid at step 500:**
```bash
python main.py --plot_step 500 --steps 2000
```

**Faster spreading disease (high sigma, lower gamma):**
```bash
python main.py --sigma 0.3 --gamma 0.002 --steps 3000
```
