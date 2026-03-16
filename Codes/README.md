# LF_MFCBF_2

MATLAB code for simulating two multi-agent scenarios in a periodic 2D domain with dangerous disks:

- `main_shepherding.m`: leader-follower / shepherding experiment with Control Barrier Functions (CBFs)
- `main_coverage.m`: diffusive coverage-style experiment with CBF safety filtering

This public repository is intended to provide the simulation code. Plotting scripts and large simulation datasets may be kept separate from the GitHub version.

## Project Structure

### Main simulation scripts

- `main_shepherding.m`
  Runs the leader-follower simulation for 50 realizations, saves trajectories and safety/performance metrics.

- `main_coverage.m`
  Runs the diffusive follower-only experiment for 50 realizations and saves trajectories and safety metrics.

### Core helper functions

- `compute_q_T.m`
  Builds the periodic interaction kernel used in the shepherding model.

- `leaders_control_macropt.m`
  Computes the leader control from the macroscopic KL-based objective.

- `compute_interactions.m`
  Computes leader-to-follower interaction velocities.

- `compute_safe_velocity_Rk.m`
  Applies the CBF-based safety filter.

- `cbf_solver.m`
  Solves the quadratic program defining the safe control action.

- `density_estimation.m`
  Estimates smooth periodic densities from particle positions.

- `generate_initial_positions.m`
  Samples admissible initial agent positions outside the dangerous disks.

- `generate_obstacle_centers.m`
  Randomly places the dangerous disk centers in the periodic domain.

- `generate_obstacles.m`
  Samples point clouds inside the dangerous disks for the kernel-based safety term.

- `gradient_2Dperiodic.m`
  Computes periodic finite-difference gradients on the square torus.

- `kernel.m`
  Backward-compatible wrapper for the Gaussian periodic kernel.

- `wrapToL.m`
  Backward-compatible wrapper for periodic wrapping into `[-L, L)`.

## Required Data

The code expects some precomputed `.mat` files to be present in the project folder, for example:

- `desired_leaders_density2.mat`

## Typical Workflow

### 1. Run simulations

For the shepherding experiment:

```matlab
main_shepherding
```

For the coverage experiment:

```matlab
main_coverage
```

The generated outputs are stored locally in the folders created by the main scripts.

## Notes

- The simulations are written for a periodic square domain.
- Dangerous regions are modeled as circular dangerous disks.
- The repository is organized around reproducible simulation runs rather than figure generation.

## MATLAB Requirements

The code uses standard MATLAB functionality and `quadprog`, so the Optimization Toolbox is required for the CBF quadratic program.
