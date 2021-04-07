import math
import numpy as np

UO2_k = 0.035            # UO2 thermal conductivity (W/cm/K)
water_cp = 4.184         # UO2 specific heat capacity (J/g/K)
htc = 0.5                # fluid-solid heat transfer coefficient (W/cm2/K)
mdot = 300.0             # mass flowrate (g/s)

# Fluid temperature model

# The fluid temperature in layer i is determined from steady-state convection
# from the solid pincell; axial conduction in the solid phase is neglected so that
# the heat deposited into the fluid layer i is the same as the heat deposited
# into the i-th solid layer. An energy balance over layer i gives:
#
# q_i=m * C_p * (T_{f,i+1/2}-T_{f,i-1/2})
#
# where 1/2 indicates temperatures on the faces of the fluid cell, m is the mass
# flowrate, and C_p is the fluid isobaric specific heat capacity. After computing
# the fluid temperatures, the fluid density is then computed as a function of this temperature.

# Solid temperature model

# The solid temperature in layer i is determined from the analytical solution
# to the heat conduction equation with constant thermal conductivity and constant heat source,
#
# -1 / r * d/dr (rk dT_{s,i}/dr )=q_i
#
# where k is the thermal conductivity. We apply boundary conditions of
# dT_{s,i}(0)/dr=0, i.e. symmetry at the pincell centerline, and
#
# q_i / (2\pi RH) = h (T_{s,i}(R)-T_{f,i})
#
# where h is a heat transfer coefficient, R is the outer radius of the pincell,
# and H is the height of the i-th layer. After solving for the solid temperature
# distribution in each layer, we simply average to obtain a single value per layer.

def water_density(T):
  """
  Returns water density at a temperature T (K) in units of g/cm3
  """
  rho = 0.001 * (0.14395 / math.pow(0.0112, 1.0 + math.pow(1.0 - T / 649.727, 0.05107)))
  return rho

def fluid_temperature(q, T_inlet, N):
  """
  Returns the fluid temperature for each layer for given solid heat source q,
  fluid inlet temperature T_inlet and number of layers N.
  """

  fluid_face_temps = np.zeros(N + 1)
  fluid_face_temps[0] = T_inlet
  for j in range(1, N + 1):
    fluid_face_temps[j] = q[j - 1] / mdot / water_cp + fluid_face_temps[j - 1]

  fluid_cell_temps = np.zeros(N)
  for j in range(N):
    fluid_cell_temps[j] = (fluid_face_temps[j] + fluid_face_temps[j + 1]) / 2.0

  return fluid_cell_temps

def fluid_density(T, N):
  """
  Returns the fluid density for each layer given a fluid temperature T and
  number of layers N.
  """

  fluid_cell_densities = np.array([water_density(i) for i in T])
  return fluid_cell_densities

def solid_temperature(q, T, N, R, Rf, H):
  """
  Returns the solid temperature for each layer given a heat source q,
  a fluid temperature T, number of layers N, pincell outer radius R,
  fuel pellet outer radius Rf, and height H.
  """
  solid_cell_temps = np.zeros(N)

  for i in range(N):
    heat_flux = q[i] / (2.0 * math.pi * R * H)
    T_solid_surface = heat_flux / htc + T[i]

    volumetric_q = q[i] / (math.pi * Rf * Rf * H)
    solid_cell_temps[i] = T_solid_surface + volumetric_q * Rf * Rf / (8 * UO2_k)

  return solid_cell_temps
