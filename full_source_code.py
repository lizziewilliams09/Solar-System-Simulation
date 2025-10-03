# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:33:42 2023

@author: eccwi
"""
from particle3D import Particle3D
import matplotlib.pyplot as pyplot
import numpy as np
import sys
from basic_functions import compute_separations
from basic_functions import compute_forces_potential 
from scipy.signal import find_peaks

def read_initial_conditions(inputfile_name):
    """

    Parameters
    ----------
    inputfile_name : string
        The name of the file with the objects and their initial conditions.

    Returns
    -------
    objects : array of Particle3D objects
        Array of the objects in the system in the form of a Particle3D .

    """
    # open file and read it
    inputfile = open(inputfile_name, "r")
    objects = []
    for line in inputfile.readlines():
        # append each Particle3D object in the array "objects"
        p = Particle3D.read_line(line)
        objects.append(p)
    return objects

def subtract_com_velocity(objects):
    """
    Parameters
    ----------
    objects : array of Particle3D objects
        Array of the objects in the system in the form of a Particle3D .

    Returns
    -------
    objects : array of Particle3D objects
        Array of the objects in the system in the form of a Particle3D, now with a centre of mass correction
        
    """
    
    #subtract the centre-of-mass velocity from a list of particles
    for j in range(len(objects)):
        objects[j].velocity = objects[j].velocity - Particle3D.com_velocity(objects)
    return objects
  
def simulation(dt, numstep, objects, outfile_name):
    """

    Parameters
    ----------
    dt : float
        time step.
    numstep : integer
        number of steps.
    objects : array of Particle3D objects
        Array of the objects in the system in the form of a Particle3D .
    outfile_name : string
        name of file to store XYZ trajectories.

    Returns
    -------
    trajectories : 3D array
        stores trajectories of each object at each time step in x, y, z directions.
    total_energy : 1D array
        total energy at each time step.

    """
    time = 0
    times = []
    #initial conditions
    separations = compute_separations(objects)
    forces, potential = compute_forces_potential(objects, separations)
    outfile = open(outfile_name, "w")
        
    #put into file
    outfile.write(f"{len(objects)}\n")
    outfile.write("Point = 0\n")
    for j in range(len(objects)):
        outfile.write(objects[j].__str__() + "\n")
    
    # Create arrays to store the results
    total_energy = np.zeros(numstep)
    trajectories = np.zeros((numstep, len(objects), 3))
    
    for i in range(numstep):

        times.append(time)
       
        #Update particle positions and update in file
        for j in range(len(objects)):
            objects[j].update_position_2nd(dt, forces[j])
        
        separations = compute_separations(objects)
        # Get the force value for the new positions
        forces_new, potential_new = compute_forces_potential(objects, separations)    
        
        # Update particle velocities by averaging current and new forces
        for j in range(len(objects)):
            objects[j].update_velocity(dt, 0.5*(forces[j]+forces_new[j]))
                 
        #put into file
        outfile.write(f"{len(objects)}\n")
        outfile.write(f"Point = {i+1} \n")
        for j in range(len(objects)):
            outfile.write(objects[j].__str__() + "\n")    
            
                            
        for j in range(len(objects)):
            forces[j] = forces_new[j] # Re-define force value for the next iteration 
            trajectories[i,j] = objects[j].position # Put locations into trajectory array  
        
        # total energy is the total potential energy and the total kinetic energy
        total_energy[i] = potential_new + Particle3D.total_kinetic_energy(objects)
        
        time += dt

    return trajectories, total_energy, times

def orbital_parameters(trajectories, dt, numstep, objects, times):
    """

    Parameters
    ----------
    trajectories : 3D array
        stores trajectories of each object at each time step in x, y, z directions.
    dt : float
        time step.
    numstep : integer
        number of steps.
    objects : array of Particle3D objects
        Array of the objects in the system in the form of a Particle3D .

    Returns
    -------
    aphelion : 1D array of floats
        Furthest distance of each object to the thing it is orbiting.
    perihelion : 1D array of floats
        Closest distance of each object to the thing it is orbiting.
    orbital_period : 1D array of floats
        time taken for object to complete a full orbit.

    """

    # Create arrays to store the results
    aphelion = np.zeros(len(objects) - 1)
    perihelion = np.zeros(len(objects) - 1)
    orbital_period = np.zeros(len(objects) - 1)
    
    # Create array to store the orbital radius of the object
    distances = np.zeros((numstep, len(objects) - 1))

    # Loop over each planet
    for j in range(len(objects) - 1):
        
        # allows for the orbit of the moon to be calculated around the earth
        # and the other objects to be calculated around the sun
        if objects[j+1].label == "Moon":
            n = 3 #earth is object 3
        else:
            n = 0 #sun is object 0
        
        # Calculate the distance of each planet/moon from the sun/earth for each time step
        for i in range(numstep):
                distances[i, j] = np.linalg.norm(trajectories[i, j+1] - trajectories[i, n])

        # Find the minimum and maximum distances
        min_distance = np.min(distances[:, j])
        max_distance = np.max(distances[:, j])
        
        # Create arrays with the positions of each peak (aphelion) and trough (perihelion)
        peaks = find_peaks(distances[:, j])[0]
        troughs = find_peaks(-distances[:, j])[0] 
        # Multiply by dt to find the times of the peaks and troughs 
        aphelion_times = np.array(peaks) * dt
        perihelion_times = np.array(troughs) * dt

        # If the planet has not completed a half orbit, return NaN for all parameters
        if len(perihelion_times) == 0 or len(aphelion_times) == 0:
            aphelion[j] = np.nan
            perihelion[j] = np.nan
            orbital_period[j] = np.nan
        # If the planets has completeled less than 1 orbit (but more than half)
        # The orbit is calculated using the difference between the perihelion and aphelion times * 2
        elif len(perihelion_times) == 1 or len(aphelion_times) == 1:
            orbital_period[j] = abs(perihelion_times[0] - aphelion_times[0])*2
            # Store the results
            perihelion[j] = min_distance
            aphelion[j] = max_distance
        else:
            # Calculate the average time between each perihelion and aphelion
            avg_perihelion_period = np.mean(np.diff(perihelion_times))
            avg_aphelion_period = np.mean(np.diff(aphelion_times))
            # The orbital period is the average of the periods above
            orbital_period[j] = (avg_perihelion_period + avg_aphelion_period)/2
            # Store the results
            perihelion[j] = min_distance
            aphelion[j] = max_distance
       
        #print information
        print(objects[j+1].label + " parameters:")
        print("Orbital Period: " + str(orbital_period[j]) + " [earth days]")     
        print("Aphelion: " + str(aphelion[j]) + "[AU]")   
        print("Perihelion: " + str(perihelion[j]) + " [AU hi]")  
        print(" ")
        
        pyplot.figure()
        pyplot.title("Distance of " + objects[j+1].label + " from " + objects[n].label)
        pyplot.xlabel("Time [earth days]")
        pyplot.ylabel("Distance [AU]")
        pyplot.plot(times, distances[:, j])
        pyplot.show()   
    
    return aphelion, perihelion, orbital_period
    
    
def main():
    dt = float(sys.argv[1]) #time step
    numstep = int(sys.argv[2]) #number of steps
    inputfile_name = sys.argv[3] #name of the input file with the initial conditions
    outfile_name = sys.argv[4] #name of output file to store XYZ trajectory
    objects_adjusted = read_initial_conditions(inputfile_name)
    trajectories, total_energy, times = simulation(dt, numstep, objects_adjusted, outfile_name)
    orbital_parameters(trajectories, dt, numstep, objects_adjusted, times)
    
if __name__ == "__main__":
    main()
    
    
  