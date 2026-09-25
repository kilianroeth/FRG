# FRG - Functional Renormalization Group Flow Solver

A numerical solver, written in C++, for functional renormalization group (FRG) flow equations. It integrates the Wetterich equation in the local potential approximation (LPA) for the effective potential of a scalar field theory, using an adaptive Runge-Kutta scheme.

Two models are currently implemented:
- **$\phi^4$ $O(N)$** theory in vacuum ```src/phi4_lpa.cpp```
- **Quark-Meson model** in vacuum ```src/quark_meson.cpp```

![effective potential flow](/results/phi4/eff_pot_flow_3.png)
*Flow of the effective potential $V_k(\rho)$ for the $\phi^4$ $O(N)$ model as the RG scale k is lowered from the UV towards the IR.*

## Background

The solver integrates the flow of the effective potential $V_k(\rho)$ with the RG "time" $t=\log(k/\Lambda)$, where $k$ is the renormalization scale and $\rho = \phi^a \phi^a /2 $ the $O(N)$ invariant field. For $\phi^4$ theory in $d=3$ dimensions this takes the form

$$
\partial_t V_k(\rho) = \frac{1}{6 \pi^2} \frac{k^3}{1 + \frac{V_k' + 2 \rho V_k''}{k^2}}
$$

The grid, spatial derivatives and RHS evaluation are model-specific; the adaptive step-doubling RK4 integrate ```src/integrator.cpp``` is shared across models.

## Build

**Prerequisites:** Cmake $\geq$3.10, Ninja, a C++20 compiler, HDF5 and OpenMP.

```
./build.sh              # debug build
./build.sh --release    # optimized build
./build.sh --run        # build and run
./build.sh --clean      # remove build/ first
```

The executable is written to ```build/``` and picked up automatically by ```--run```.

## Results

TODO

## Plotting

Plot and post-processing live in ```plot/(phi4.ipynb,QM.ipynb)```. Install the Python dependencies with:

```
pip install -r requirements.txt
```

## Repository layout

|direcotry|content|
|---|---|
|stc|solver implementation (grid, integrator, physics models, main)|
|include| headers|
|plot| Python plotting notebooks/scripts|
|results|outut data and figures (mostly gitignored)|