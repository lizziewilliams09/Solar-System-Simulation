# Solar System Simulation with Velocity Verlet Integration

**Description:**
This project simulates the motion of planets, the Moon, and comets in the Solar System using the Velocity Verlet integration method in Python. The simulation calculates trajectories, total energy, and orbital parameters such as aphelion, perihelion, and orbital period.

The analysis of observables (orbital parameters, energy conservation, convergence) was performed in Excel and is included in the accompanying report.

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

* The Excel analysis is not included in the repo; it is summarized in the report.
* The simulation assumes point-mass particles and does not include relativistic corrections.
