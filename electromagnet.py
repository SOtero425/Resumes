"""
Created on Sat Jun 18 19:50:11 2022

@author: sotero425
"""
import scipy as sp
import numpy as np

'''
Library Outline
-Class for magnetic field (B)
    -Attributes: Shape/Dimensions, B-field at each point,
    interpolation function?

-Class for electromagnets
    -Attributes: Layers, turns/layer, wire thickness.
        -Initially take single values, but eventually add ability to input
        arrays.

-Method to calculate magnetic field from electromagnet at a specified point
    in space (Biot-savart law), probably based off current = 1 amp (linear
    scaling can be applied for different current values afterwards).

-Possibly a method to generate the magnetic field from electromagnet within a
    2D grid.
    -Interpolate at this point or wait until afterwards? Is it more efficient,
    computationally, to interpolate at this point or after summing magnetic
    fields together?

-Possibly the interpolation function.

-Function to expand the single-quadrant grid values to 2-quadrant.

-Function to calculate field at any point in a cylinder utilizing cylindrical
    symmetry.

-Function to sum fields from multiple coils at different points in space.

-F=q*vxB?

-Runge-kutta for velocity, and position calculations?

-Class for charged particle?
    -Attributes: Charge, mass, velocity, position.
'''

'''
Two functions: One to calculate the force on a charged particle with scalar and
angular information, a second that requires velocity and position in cartesian
vector form.
'''


def force_simple(charge=1,
                 velocity=1,
                 bfield=1,
                 theta=0):

    force = charge * velocity * bfield * np.cos()
    return force


def force(charge=1,
          velocity=np.array([1, 1, 1]),
          bfield=np.array([1, 1, 1])):

    force = charge * np.cross(velocity, bfield)
    return force


class BField():
    def __init__(self,
                 shape=1):

        self.shape = np.zeros([shape])


class Electromagnet():
    def __init__(self,
                 layers=1,
                 turns=1,
                 wire_thickness=1,
                 current=1,
                 radius=1):

        self.layers = np.int(layers)
        self.turns = np.int(turns)
        self.wire_thickness = np.int(wire_thickness)
        self.current = np.int(current)
        self.radius = np.int(radius)
        self.__constant = (sp.m_0 * self.current / (4*sp.pi))


# Biot-savart
# dB = (m_0*I/4*pi*r^2)*(dL x r')/r'^3
# Integrate dB using Simpson's Rule then return B
# Need to convert dL segments to cartesian coordinates

    def bfield(self,
               position=np.array([0, 0, 0]),
               slices=np.int(99)):

        position = np.array(position)
        slices = np.int(slices)
        d_theta = 2*sp.pi/slices  # in radians
        theta = 0
        dB = np.zeros([slices, 3])

        for ind in range(slices):
            x, y = self.radius*np.cos(theta), self.radius*np.sin(theta)
            dL = [x, y, 0]
            theta += d_theta
            r_prime = position - dL
            r_prime_mag = np.linalg.norm(r_prime)
            dB[ind] = self.__constant * np.cross(dL, r_prime) / r_prime_mag**3

        B = sp.integrate.simpson(dB, x=None, dx=d_theta)
        return B


class ChargedParticle():
    def __init__(self,
                 charge=1,
                 mass=1,
                 velocity=1,
                 position=1):

        self.charge = charge
        self.mass = mass
        self.velocity = velocity
        self.position = position
