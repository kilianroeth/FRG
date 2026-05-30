# FRG - Functional Renormalization Group

Trying to implement some FRG numerics from scratch.

Numerical implementation of the Wetterich equation in the Loacl Potential Approximation for $\phi^4$ theory in $d=3$ dimensions.

## Physics

We want to solve the flow equation for the effective potential $V_k(\rho)$

$$
\partial_t V_k(\rho) = \frac{1}{6 \pi^2} \frac{k^3}{1+ \frac{V'_k(\rho) + 2 \rho V''_k(\rho)}{k^2}}
$$

where $t = \log(k/\Lambda)$ is the RG-time and $\rho = \frac{1}{2} \phi^2$ is the $O(N)$ invariant field. The flow runs from the cutoff scale $\Lambda$ where $V_{\Lambda}$ is the classical potential and all quantum fluctuations are suppressed to $k=0$ where all quantum fluctutations are present. (UV: t = 0, IR: t $\to -\infty$).

The classcial initial contition at $k = \Lambda$ is

$$
V_\Lambda(\rho) = m_\Lambda^2 \rho + \frac{\lambda_\Lambda}{3!} \rho^2
$$

## Parameters

| Paramter | |
| ---------| - |
| `m2` | Mass paramter $m^2$ at the UV scale. Negative $\to$ broken phase |
| `lamda` | Quartic coupling at the UV scale. |
| `t_start` | RG time in the UV (t_start = 0) |
| `t_end` | RG time in the IR (eg. -15.0 so $k_{IR} = e^{-15}$) |
| `n_rho` | Number of grid points in rho |
| `rho_max` | Upper boundary of the rho grid (lower bound is always 0) |

## build & Run

```bash
./build.sh              # Build (Debug) 
./build.sh --run        # Build and run
./build.sh --release    # Build in Release mode
./build.sh --clean      # Clean build
```
## Output

Results are written to `results/.` The main output `flow(_adaptive).csv` contains two blocks: `# block: V` with $V_k(\rho)$ at each snapshot $k$, and `# block: RHS` with the corresponding flow equation RHS. Use `plot/plot.ipynb` to visualize.
