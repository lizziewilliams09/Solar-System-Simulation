# -*- coding: utf-8 -*-
"""
First Submission: Basic Functions
"""

import numpy as np
from particle3D import Particle3D

def compute_separations(particles):
    """
    Parameters
    ----------
    particles : list
        Particles in the system.

    Returns
    -------
    separations : 3d array
        Vector deparation of each pair of particles.

    """
    n = len(particles)
    separations = np.zeros((n, n, 3)) #create array to store the separations
    for i in range(0, n): #goes through each particle
        position_i = particles[i].position
        for j in range (i+1,n): #find separation with each other particle without double counting
            position_j = particles[j].position
            r = np.subtract(position_i, position_j) #find seperation vector r
            #add seperations to array
            separations[i,j] = r
            separations[j,i] = -r #the vector between i and j is the negative of the vector between j and i
    return separations
    
def compute_forces_potential(particles, separations):
    """
    Parameters
    ----------
    particles : list
        Particles in the system.
    separations : 3d array
        Vector deparation of each pair of particles.

    Returns
    -------
    forces : 2d array
        Force on each particle.
    potential : float
        Total system potential energy.

    """
    n = len(particles)
    G = 8.887692593*10**(-10) #gravitational constant for force equation
    forces = np.zeros((n,3)) #create an array to store forces on all particles
    #forces
    for i in range(0,n): #run through each particle
        force_i = np.zeros(3) #create an array to store the force on particle i
        for j in range (0,n): #find the force between i and the other particles, j
            r = separations[i,j] #gives seperation vector r between i and j
            modulus_r = np.linalg.norm(r) #modulus of r
            if (i != j): 
                #add force between i and j to the array of forces on i
                force_i += (-G*particles[i].mass*particles[j].mass*r)/(modulus_r**3) 
        forces[i] = force_i #add forces on i to complete forces array
    #potentials
    potential = 0 #create a float to store the potential
    for i in range(0, n): #run through each particle
        for j in range (i+1,n): #run through each particle after particle i (to avoid double counting)
            modulus_r = np.linalg.norm([separations[i,j]]) #gives modulus of seperation vector r between i and j
            potential += (-G*particles[i].mass*particles[j].mass)/(modulus_r) #add potential between i and j to potential
    return forces, potential
