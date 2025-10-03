# Solar System Simulation with Velocity Verlet Integration

**Description:**
Simulation of the Solar System using the Velocity Verlet integration scheme in Python, modelling planets, the Moon, Pluto, and Halley’s Comet as gravitational N-body point particles. The code outputs orbital parameters (perihelion, aphelion, orbital periods) and tracks energy conservation. 

Convergence testing identified an optimal time-step (0.5 days) for ≤0.5% relative error. Observables were validated against real data, including the Moon and Jupiter, with most errors <2%, and Neptune’s anomaly analyzed. Kepler’s Third Law was verified, including a “super-Jupiter” simulation to illustrate mass effects on orbital dynamics. All trajectory and observable analyses, convergence tests, and verification steps were performed in Excel and are fully documented in the accompanying report.

---

## **Files and Structure**

* `particle3D.py` - Class defining 3D particles, including position, velocity, kinetic energy, and update methods.
* `basic_functions.py` - Functions to compute inter-particle separations, forces, and potential energy.
* `full_source_code.py` - Main simulation script. Reads initial conditions, runs simulation, computes orbital parameters, and plots distances.
* `solar_system.txt` - Initial conditions for the Sun, planets, Moon, Pluto, and Halley’s comet.
* `XYZ_file.txt` - Sample output of object trajectories at each time step.
* `Simulating the Solar System Report.pdf` - PDF report describing methodology, results, and analysis.


---

## **Usage**

1. Ensure all `.py` files and the `solar_system.txt` initial conditions file are in the same directory.
2. Run the simulation from the command line:

```bash
python full_source_code.py <dt> <numsteps> solar_system.txt XYZ_file.txt
```

* `<dt>` = timestep (e.g., 0.01)
* `<numsteps>` = number of simulation steps (e.g., 1000)

3. Trajectories will be saved in `XYZ_file.txt`.
4. Orbital parameters are printed to the console, and distance plots are displayed interactively.
5. See **Simulating the Solar System Report.pdf** for full analysis, plots, and discussion.

---

## **Notes**

* The Excel analysis is not included in the repo; it is summarised in the report.
* The simulation assumes point-mass particles and does not include relativistic corrections.
