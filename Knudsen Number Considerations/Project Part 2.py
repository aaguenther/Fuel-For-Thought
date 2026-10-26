import math as m
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Given

T_initial   = 298 # K
P_c_bar     = np.array([1, 10, 100]) # bar
P_c         = P_c_bar * 1e5 # Pa
phi         = [3.0]
R_u         = ct.gas_constant # J / (kmol * K)
N_A         = 6.0221408e23
k_b         = 1.38064852e-23
theta       = 15 * (m.pi / 180) # radians
Pa_to_atm   = 9.86923e-6        # atm/Pa


# Have numpy print out numbers using scientific notation
np.set_printoptions(suppress=False, precision=4)


# Calculate equilibrium state in the combustor
gas_1bar = ct.Solution("gri30.yaml") # Mechanism for H2
gas_10bar = ct.Solution("gri30.yaml")
gas_100bar = ct.Solution("gri30.yaml")

hydrogen = "H2"
oxygen = "O2"

from src.equilibrium import equilibrium_properties

chamber_properties_1_bar = equilibrium_properties(hydrogen, oxygen, T_initial, P_c[0], phi, gas_1bar)
chamber_properties_10_bar = equilibrium_properties(hydrogen, oxygen, T_initial, P_c[1], phi, gas_10bar)
chamber_properties_100_bar = equilibrium_properties(hydrogen, oxygen, T_initial, P_c[2], phi, gas_100bar)

# Store equilibrium state properties in arrays

T_c = [chamber_properties_1_bar[0][0], chamber_properties_10_bar[0][0], chamber_properties_100_bar[0][0]]
gamma = [chamber_properties_1_bar[1][0], chamber_properties_10_bar[1][0], chamber_properties_100_bar[1][0]]
R = [R_u / chamber_properties_1_bar[2][0], R_u / chamber_properties_10_bar[2][0], R_u / chamber_properties_100_bar[2][0]]

# Mean molecular weight is in units kg/molecule
mw = [chamber_properties_1_bar[2][0] / (1000 * N_A), chamber_properties_10_bar[2][0] / (1000* N_A), chamber_properties_100_bar[2][0] / (1000 * N_A)]

# Calculate c* for each chamber pressure

from src.equilibrium import characteristic_velocity

c_star_1bar = characteristic_velocity(T_c[0], gamma[0], R[0])
c_star_10bar = characteristic_velocity(T_c[1], gamma[1], R[1])
c_star_100bar = characteristic_velocity(T_c[2], gamma[2], R[2])

# Calculate mu, d, and lambda for each chamber pressure

mu_1_bar_chamber = gas_1bar.viscosity
mu_10_bar_chamber = gas_10bar.viscosity
mu_100_bar_chamber = gas_100bar.viscosity

d_1_bar_chamber = m.sqrt((5 / 16) * (m.sqrt(m.pi * mw[0] * k_b * T_c[0]) / (m.pi * mu_1_bar_chamber)))
d_10_bar_chamber = m.sqrt((5 / 16) * (m.sqrt(m.pi * mw[1] * k_b * T_c[1]) / (m.pi * mu_10_bar_chamber)))
d_100_bar_chamber = m.sqrt((5 / 16) * (m.sqrt(m.pi * mw[2] * k_b * T_c[2]) / (m.pi * mu_100_bar_chamber)))

n_1_bar_chamber = P_c[0] / (k_b * T_c[0])
n_10_bar_chamber = P_c[1] / (k_b * T_c[1])
n_100_bar_chamber = P_c[2] / (k_b * T_c[2])

lambda_1_bar_chamber = 1 / (m.sqrt(2) * n_1_bar_chamber * m.pi * d_1_bar_chamber ** 2)
lambda_10_bar_chamber = 1 / (m.sqrt(2) * n_10_bar_chamber * m.pi * d_10_bar_chamber ** 2)
lambda_100_bar_chamber = 1 / (m.sqrt(2) * n_100_bar_chamber * m.pi * d_100_bar_chamber ** 2)

############################################# TASK 1 ###################################################################

# Calculate properties along diverging section

# Create arrays

area_ratio = np.linspace(1, 50, 1000)


M_1_bar = []
M_10_bar = []
M_100_bar = []

rho_1_bar = []
rho_10_bar = []
rho_100_bar = []

T_1_bar = []
T_10_bar = []
T_100_bar = []

P_1_bar = []
P_10_bar = []
P_100_bar = []

u_1_bar = []
u_10_bar = []
u_100_bar = []

mu_1_bar = []
mu_10_bar = []
mu_100_bar = []

d_1_bar = []
d_10_bar = []
d_100_bar = []

n_1_bar = []
n_10_bar = []
n_100_bar = []

lambda_1_bar = []
lambda_10_bar = []
lambda_100_bar = []

# Populate arrays

from src.quasi_1D_flows import area_mach_relation
from src.isentropic_relations import po_p
from src.isentropic_relations import to_t

for j in area_ratio:

    M_1_bar_j = area_mach_relation(j, gamma[0], "supersonic")
    M_10_bar_j = area_mach_relation(j, gamma[1], "supersonic")
    M_100_bar_j = area_mach_relation(j, gamma[2], "supersonic")

    M_1_bar.append(M_1_bar_j)
    M_10_bar.append(M_10_bar_j)
    M_100_bar.append(M_100_bar_j)

    P_1_bar_j = (1 / po_p(gamma[0], M_1_bar_j)) * (P_c[0])
    P_10_bar_j = (1 / po_p(gamma[1], M_10_bar_j)) * (P_c[1])
    P_100_bar_j = (1 / po_p(gamma[2], M_100_bar_j)) * (P_c[2])

    P_1_bar.append(P_1_bar_j)
    P_10_bar.append(P_10_bar_j)
    P_100_bar.append(P_100_bar_j)

    T_1_bar_j = (1 / to_t(gamma[0], M_1_bar_j)) * T_c[0]
    T_10_bar_j = (1 / to_t(gamma[1], M_10_bar_j)) * T_c[1]
    T_100_bar_j = (1 / to_t(gamma[2], M_100_bar_j)) * T_c[2]

    T_1_bar.append(T_1_bar_j)
    T_10_bar.append(T_10_bar_j)
    T_100_bar.append(T_100_bar_j)

    rho_1_bar_j = P_1_bar_j / (R[0] * T_1_bar_j)
    rho_10_bar_j = P_10_bar_j / (R[1] * T_10_bar_j)
    rho_100_bar_j = P_100_bar_j / (R[2] * T_100_bar_j)

    rho_1_bar.append(rho_1_bar_j)
    rho_10_bar.append(rho_10_bar_j)
    rho_100_bar.append(rho_100_bar_j)

    u_1_bar_j = M_1_bar_j * m.sqrt(gamma[0] * R[0] * T_1_bar_j)
    u_10_bar_j = M_10_bar_j * m.sqrt(gamma[1] * R[1] * T_10_bar_j)
    u_100_bar_j = M_100_bar_j * m.sqrt(gamma[2] * R[2] * T_100_bar_j)

    u_1_bar.append(u_1_bar_j)
    u_10_bar.append(u_10_bar_j)
    u_100_bar.append(u_100_bar_j)

    # Set gas objects containing combustion products to the local temp and pressure in the nozzle

    gas_1bar.TP = T_1_bar_j, P_1_bar_j
    gas_10bar.TP = T_10_bar_j, P_10_bar_j
    gas_100bar.TP = T_100_bar_j, P_100_bar_j

    # Calculate viscosity

    mu_1_bar_j = gas_1bar.viscosity
    mu_10_bar_j = gas_10bar.viscosity
    mu_100_bar_j = gas_100bar.viscosity

    mu_1_bar.append(mu_1_bar_j)
    mu_10_bar.append(mu_10_bar_j)
    mu_100_bar.append(mu_100_bar_j)

    # Calculate effective diameter

    d_1_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[0] * k_b * T_1_bar_j) / (m.pi * mu_1_bar_j)))
    d_10_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[1] * k_b * T_10_bar_j) / (m.pi * mu_10_bar_j)))
    d_100_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[2] * k_b * T_100_bar_j) / (m.pi * mu_100_bar_j)))

    d_1_bar.append(d_1_bar_j)
    d_10_bar.append(d_10_bar_j)
    d_100_bar.append(d_100_bar_j)

    # Calculate number density

    n_1_bar_j = P_1_bar_j / (k_b * T_1_bar_j)
    n_10_bar_j = P_10_bar_j / (k_b * T_10_bar_j)
    n_100_bar_j = P_100_bar_j / (k_b * T_100_bar_j)

    n_1_bar.append(n_1_bar_j)
    n_10_bar.append(n_10_bar_j)
    n_100_bar.append(n_100_bar_j)

    lambda_1_bar_j = 1 / (m.sqrt(2) * n_1_bar_j * m.pi * d_1_bar_j**2)
    lambda_10_bar_j = 1 / (m.sqrt(2) * n_10_bar_j * m.pi * d_10_bar_j**2)
    lambda_100_bar_j = 1 / (m.sqrt(2) * n_100_bar_j * m.pi * d_100_bar_j**2)

    lambda_1_bar.append(lambda_1_bar_j)
    lambda_10_bar.append(lambda_10_bar_j)
    lambda_100_bar.append(lambda_100_bar_j)

############################################# TASK 2 ###################################################################

# Given

C_F = 1.9

# Create array of Knudsen numbers

Kn = np.linspace(.001, 1, 1000)

# Create empty arrays for throat and exit areas

A_t_1bar = []
A_t_10bar = []
A_t_100bar = []
A_e_1bar = []
A_e_10bar = []
A_e_100bar = []


for i in Kn:

    # Calculate throat diameter to satisfy each Kn
    D_t_1bar_i = lambda_1_bar[0] / i
    D_t_10bar_i = lambda_10_bar[0] / i
    D_t_100bar_i = lambda_100_bar[0] / i

    # Calculate corresponding area
    A_t_1bar_i = 0.25 * m.pi * D_t_1bar_i**2
    A_t_10bar_i = 0.25 * m.pi * D_t_10bar_i**2
    A_t_100bar_i = 0.25 * m.pi * D_t_100bar_i**2

    A_t_1bar.append(A_t_1bar_i)
    A_t_10bar.append(A_t_10bar_i)
    A_t_100bar.append(A_t_100bar_i)

    # Calculate exit diameter to satisfy each Kn
    D_e_1bar_i = lambda_1_bar[len(area_ratio) - 1] / i
    D_e_10bar_i = lambda_10_bar[len(area_ratio) - 1] / i
    D_e_100bar_i = lambda_100_bar[len(area_ratio) - 1] / i

    # Calculate corresponding area
    A_e_1bar_i = 0.25 * m.pi * D_e_1bar_i**2
    A_e_10bar_i = 0.25 * m.pi * D_e_10bar_i**2
    A_e_100bar_i = 0.25 * m.pi * D_e_100bar_i**2

    A_e_1bar.append(A_e_1bar_i)
    A_e_10bar.append(A_e_10bar_i)
    A_e_100bar.append(A_e_100bar_i)

# Select a design thrust

T = 50 # mN
T_N = T/1000 # Convert to newtons

# Find design throat and exit areas

A_t_design_1bar = T_N / (C_F * P_c[0])
A_t_design_10bar = T_N / (C_F * P_c[1])
A_t_design_100bar = T_N / (C_F * P_c[2])

A_e_design_1bar = 50 * A_t_design_1bar
A_e_design_10bar = 50 * A_t_design_10bar
A_e_design_100bar = 50 * A_t_design_100bar

# Find corresponding Kn values

Kn_t_design_1bar = lambda_1_bar[0] / (m.sqrt((4 * A_t_design_1bar) / m.pi))
Kn_t_design_10bar = lambda_10_bar[0] / (m.sqrt((4 * A_t_design_10bar) / m.pi))
Kn_t_design_100bar = lambda_100_bar[0] / (m.sqrt((4 * A_t_design_100bar) / m.pi))

Kn_e_design_1bar = lambda_1_bar[len(area_ratio) - 1] / (m.sqrt((4 * A_e_design_1bar) / m.pi))
Kn_e_design_10bar = lambda_10_bar[len(area_ratio) - 1] / (m.sqrt((4 * A_e_design_10bar) / m.pi))
Kn_e_design_100bar = lambda_100_bar[len(area_ratio) - 1] / (m.sqrt((4 * A_e_design_100bar) / m.pi))

# Find minimum thrust values for continuum

T_min_1bar = 0.25 * m.pi * (lambda_1_bar[len(area_ratio) - 1] / 0.1)**2 * C_F * P_c[0]
T_min_10bar = 0.25 * m.pi * (lambda_10_bar[len(area_ratio) - 1] / 0.1)**2 * C_F * P_c[1]
T_min_100bar = 0.25 * m.pi * (lambda_100_bar[len(area_ratio) - 1] / 0.1)**2 * C_F * P_c[2]

# Create plot
plt.plot(Kn, A_t_1bar, label = "Throat Area Required @ 1 bar", color = "red")
plt.plot(Kn, A_t_10bar, label = "Throat Area Required @ 10 bar", color = "red", linestyle = "--")
plt.plot(Kn, A_t_100bar, label = "Throat Area Required @ 100 bar", color = "red", linestyle = "-.")
plt.plot(Kn, A_e_1bar, label = "Exit Area Required @ 1 bar", color = "blue")
plt.plot(Kn, A_e_10bar, label = "Exit Area Required @ 10 bar", color = "blue", linestyle = "--")
plt.plot(Kn, A_e_100bar, label = "Exit Area Required @ 100 bar", color = "blue", linestyle = "-.")

# Add design points

plt.scatter(Kn_t_design_1bar, A_t_design_1bar, color = "red")
plt.scatter(Kn_t_design_10bar, A_t_design_10bar, color = "red")
plt.scatter(Kn_t_design_100bar, A_t_design_100bar, color = "red")

plt.scatter(Kn_e_design_1bar, A_e_design_1bar, color = "blue")
plt.scatter(Kn_e_design_10bar, A_e_design_10bar, color = "blue")
plt.scatter(Kn_e_design_100bar, A_e_design_100bar, color = "blue")

# Format plot

plt.title("Nozzle Area Required vs. Knudsen Number")
plt.xlabel("Knudsen Number (Kn)")
plt.ylabel("Area Required (m^2)")
plt.xscale("log")
plt.yscale("log")
plt.axvspan(0, 0.01, color="grey", alpha=0.2, label="No Slip Continuum")
plt.axvspan(0.01, 0.1, color="grey", alpha=0.4, label="Slip Flow")
plt.axvspan(0.1, 1.0, color="grey", alpha=0.6, label="Transition Flow")

# Display plot

plt.legend()
plt.show()

# Print Results


print("################################ Design Parameters ####################################\n")
print("Thurst =", T, "mN\n")
print("################################ Chamber Properties ###################################\n")
print("**1 bar**\n")
print("gamma =", f"{gamma[0]:.3f}")
print("R =", f"{R[0]:.2f}", "J/(kg*K)")
print("T =", f"{T_c[0]:.2f}", "K")
print("P =", f"{P_c_bar[0]:.2f}", "bar")
print("mu =", f"{mu_1_bar_chamber:.6e}", "Pa*s")
print("d =", f"{d_1_bar_chamber:.6e}", "m")
print("lambda =", f"{lambda_1_bar_chamber:.6e}", "m")
print("c* =", f"{c_star_1bar:.2f}", "m/s\n")
print("**10 bar**\n")
print("gamma =", f"{gamma[1]:.3f}")
print("R =", f"{R[1]:.2f}", "J/(kg*K)")
print("T =", f"{T_c[1]:.2f}", "K")
print("P =", f"{P_c_bar[1]:.2f}", "bar")
print("mu =", f"{mu_10_bar_chamber:.6e}", "Pa*s")
print("d =", f"{d_10_bar_chamber:.6e}", "m")
print("lambda =", f"{lambda_10_bar_chamber:.6e}", "m")
print("c* =", f"{c_star_10bar:.2f}", "m/s\n")
print("**100 bar**\n")
print("gamma =", f"{gamma[2]:.3f}")
print("R =", f"{R[2]:.2f}", "J/(kg*K)")
print("T =", f"{T_c[2]:.2f}", "K")
print("P =", f"{P_c_bar[2]:.2f}", "bar")
print("mu =", f"{mu_100_bar_chamber:.6e}", "Pa*s")
print("d =", f"{d_100_bar_chamber:.6e}", "m")
print("lambda =", f"{lambda_100_bar_chamber:.6e}", "m")
print("c* =", f"{c_star_100bar:.2f}", "m/s\n")
print("################################ Throat Properties ###################################\n")
print("**1 bar**\n")
print("T =", f"{T_1_bar[0]:.2f}", "K")
print("P =", f"{P_1_bar[0]/1e5:.6f}", "bar")
print("mu =", f"{mu_1_bar[0]:.6e}", "Pa*s")
print("d =", f"{d_1_bar[0]:.6e}", "m")
print("lambda =", f"{lambda_1_bar[0]:.6e}", "m\n")
print("**10 bar**\n")
print("T =", f"{T_10_bar[0]:.2f}", "K")
print("P =", f"{P_10_bar[0]/1e5:.6f}", "bar")
print("mu =", f"{mu_10_bar[0]:.6e}", "Pa*s")
print("d =", f"{d_10_bar[0]:.6e}", "m")
print("lambda =", f"{lambda_10_bar[0]:.6e}", "m\n")
print("**100 bar**\n")
print("T =", f"{T_100_bar[0]:.2f}", "K")
print("P =", f"{P_100_bar[0]/1e5:.6f}", "bar")
print("mu =", f"{mu_100_bar[0]:.6e}", "Pa*s")
print("d =", f"{d_100_bar[0]:.6e}", "m")
print("lambda =", f"{lambda_100_bar[0]:.6e}", "m\n")
print("################################ Exit Properties ###################################\n")
print("**1 bar**\n")
print("T =", f"{T_1_bar[len(area_ratio) - 1]:.2f}", "K")
print("P =", f"{P_1_bar[len(area_ratio) - 1]/1e5:.6f}", "bar")
print("mu =", f"{mu_1_bar[len(area_ratio) - 1]:.6e}", "Pa*s")
print("d =", f"{d_1_bar[len(area_ratio) - 1]:.6e}", "m")
print("lambda =", f"{lambda_1_bar[len(area_ratio) - 1]:.6e}", "m\n")
print("**10 bar**\n")
print("T =", f"{T_10_bar[len(area_ratio) - 1]:.2f}", "K")
print("P =", f"{P_10_bar[len(area_ratio) - 1]/1e5:.6f}", "bar")
print("mu =", f"{mu_10_bar[len(area_ratio) - 1]:.6e}", "Pa*s")
print("d =", f"{d_10_bar[len(area_ratio) - 1]:.6e}", "m")
print("lambda =", f"{lambda_10_bar[len(area_ratio) - 1]:.6e}", "m\n")
print("**100 bar**\n")
print("T =", f"{T_100_bar[len(area_ratio) - 1]:.2f}", "K")
print("P =", f"{P_100_bar[len(area_ratio) - 1]/1e5:.6f}", "bar")
print("mu =", f"{mu_100_bar[len(area_ratio) - 1]:.6e}", "Pa*s")
print("d =", f"{d_100_bar[len(area_ratio) - 1]:.6e}", "m")
print("lambda =", f"{lambda_100_bar[len(area_ratio) - 1]:.6e}", "m\n")
print("################################ Minimum Thurst for Continuum ###################################\n")
print("T_min at 1 bar =", T_min_1bar * 1000, "mN")
print("T_min at 10 bar =", T_min_10bar * 1000, "mN")
print("T_min at 100 bar =", T_min_100bar * 1000, "mN")


print("\n")
print("################################### Task 3 Related Outputs ######################################\n")

############################################# TASK 3 ###################################################################


# Part (a) - Plot the local Knudsen nu,ber along the diverging section of the nozzle.


# Import functions from the src folder
from src.task_3_functs import x_to_Area_Ratio
from src.task_3_functs import nozzle_length
from src.task_3_functs import x_to_Diameter
from src.task_3_functs import Area_ratio_to_x
from src.task_3_functs import exp_time_funct
from src.task_3_functs import relax_time_vib


# Calculate properties along diverging section using a slightly modified area ratio array
# First, we must create the array of the x-values along the diverging nozzle


# Create an array storing the design diameters for each chamber pressure
At_vals = np.array([A_t_design_1bar,A_t_design_10bar,A_t_design_100bar])
Ae_vals = np.array([A_e_design_1bar,A_e_design_10bar,A_e_design_100bar])


# Obtain the length of the diverging nozzle section for each chamber pressure.
L = []
L.append(nozzle_length(Ae_vals[0],At_vals[0],theta))
L.append(nozzle_length(Ae_vals[1],At_vals[1],theta))
L.append(nozzle_length(Ae_vals[2],At_vals[2],theta))


# Create an array of x values linearly spanning the diverging nozzle for each pressure
x_values_1bar   = np.linspace(0, L[0], 1000)
x_values_10bar  = np.linspace(0, L[1], 1000)
x_values_100bar = np.linspace(0, L[2], 1000)


# Create an array of non-dimensionalized x-values for plotting x/Dt
x_vals_nondimen_1_bar = x_values_1bar / (np.sqrt(4*At_vals[0]/m.pi))
x_vals_nondimen_10_bar = x_values_10bar / (np.sqrt(4*At_vals[1]/m.pi))
x_vals_nondimen_100_bar = x_values_100bar / (np.sqrt(4*At_vals[2]/m.pi))


# Store the delta x values for each chamber pressure
del_x_1_bar   = x_values_1bar[1]
del_x_10_bar  = x_values_10bar[1]
del_x_100_bar = x_values_100bar[1]


# Convert these x values arrays into their corresponding area ratio arrays
Area_ratios_1bar   = x_to_Area_Ratio(At_vals[0],x_values_1bar,theta)
Area_ratios_10bar  = x_to_Area_Ratio(At_vals[1],x_values_10bar,theta)
Area_ratios_100bar = x_to_Area_Ratio(At_vals[2],x_values_100bar,theta)


# Create arrays

M_1_bar = []
M_10_bar = []
M_100_bar = []

rho_1_bar = []
rho_10_bar = []
rho_100_bar = []

T_1_bar = []
T_10_bar = []
T_100_bar = []

P_1_bar = []
P_10_bar = []
P_100_bar = []

u_1_bar = []
u_10_bar = []
u_100_bar = []

mu_1_bar = []
mu_10_bar = []
mu_100_bar = []

d_1_bar = []
d_10_bar = []
d_100_bar = []

n_1_bar = []
n_10_bar = []
n_100_bar = []

lambda_1_bar = []
lambda_10_bar = []
lambda_100_bar = []

knudsen_local_1_bar = []
knudsen_local_10_bar = []
knudsen_local_100_bar = []

nozzle_dimater_1_bar   = []
nozzle_dimater_10_bar  = []
nozzle_dimater_100_bar = []

mean_therm_speed_1_bar   = []
mean_therm_speed_10_bar  = []
mean_therm_speed_100_bar = []

collision_rate_local_1_bar   = []
collision_rate_local_10_bar  = []
collision_rate_local_100_bar = []

cummulative_collisions_1_bar   = [0]
cummulative_collisions_10_bar  = [0]
cummulative_collisions_100_bar = [0]


from src.quasi_1D_flows import area_mach_relation
from src.isentropic_relations import po_p
from src.isentropic_relations import to_t


# Use for loop to populate important properties
for j in range(0,len(Area_ratios_1bar)):

    M_1_bar_j = area_mach_relation(Area_ratios_1bar[j], gamma[0], "supersonic")
    M_10_bar_j = area_mach_relation(Area_ratios_10bar[j], gamma[1], "supersonic")
    M_100_bar_j = area_mach_relation(Area_ratios_100bar[j], gamma[2], "supersonic")

    M_1_bar.append(M_1_bar_j)
    M_10_bar.append(M_10_bar_j)
    M_100_bar.append(M_100_bar_j)

    P_1_bar_j = (1 / po_p(gamma[0], M_1_bar_j)) * (P_c[0])
    P_10_bar_j = (1 / po_p(gamma[1], M_10_bar_j)) * (P_c[1])
    P_100_bar_j = (1 / po_p(gamma[2], M_100_bar_j)) * (P_c[2])

    P_1_bar.append(P_1_bar_j)
    P_10_bar.append(P_10_bar_j)
    P_100_bar.append(P_100_bar_j)

    T_1_bar_j = (1 / to_t(gamma[0], M_1_bar_j)) * T_c[0]
    T_10_bar_j = (1 / to_t(gamma[1], M_10_bar_j)) * T_c[1]
    T_100_bar_j = (1 / to_t(gamma[2], M_100_bar_j)) * T_c[2]

    T_1_bar.append(T_1_bar_j)
    T_10_bar.append(T_10_bar_j)
    T_100_bar.append(T_100_bar_j)

    rho_1_bar_j = P_1_bar_j / (R[0] * T_1_bar_j)
    rho_10_bar_j = P_10_bar_j / (R[1] * T_10_bar_j)
    rho_100_bar_j = P_100_bar_j / (R[2] * T_100_bar_j)

    rho_1_bar.append(rho_1_bar_j)
    rho_10_bar.append(rho_10_bar_j)
    rho_100_bar.append(rho_100_bar_j)

    u_1_bar_j = M_1_bar_j * m.sqrt(gamma[0] * R[0] * T_1_bar_j)
    u_10_bar_j = M_10_bar_j * m.sqrt(gamma[1] * R[1] * T_10_bar_j)
    u_100_bar_j = M_100_bar_j * m.sqrt(gamma[2] * R[2] * T_100_bar_j)

    u_1_bar.append(u_1_bar_j)
    u_10_bar.append(u_10_bar_j)
    u_100_bar.append(u_100_bar_j)

    # Set gas objects containing combustion products to the local temp and pressure in the nozzle

    gas_1bar.TP = T_1_bar_j, P_1_bar_j
    gas_10bar.TP = T_10_bar_j, P_10_bar_j
    gas_100bar.TP = T_100_bar_j, P_100_bar_j

    # Calculate viscosity

    mu_1_bar_j = gas_1bar.viscosity
    mu_10_bar_j = gas_10bar.viscosity
    mu_100_bar_j = gas_100bar.viscosity

    mu_1_bar.append(mu_1_bar_j)
    mu_10_bar.append(mu_10_bar_j)
    mu_100_bar.append(mu_100_bar_j)

    # Calculate effective diameter

    d_1_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[0] * k_b * T_1_bar_j) / (m.pi * mu_1_bar_j)))
    d_10_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[1] * k_b * T_10_bar_j) / (m.pi * mu_10_bar_j)))
    d_100_bar_j = m.sqrt((5/16) * (m.sqrt(m.pi * mw[2] * k_b * T_100_bar_j) / (m.pi * mu_100_bar_j)))

    d_1_bar.append(d_1_bar_j)
    d_10_bar.append(d_10_bar_j)
    d_100_bar.append(d_100_bar_j)

    # Calculate number density

    n_1_bar_j = P_1_bar_j / (k_b * T_1_bar_j)
    n_10_bar_j = P_10_bar_j / (k_b * T_10_bar_j)
    n_100_bar_j = P_100_bar_j / (k_b * T_100_bar_j)

    n_1_bar.append(n_1_bar_j)
    n_10_bar.append(n_10_bar_j)
    n_100_bar.append(n_100_bar_j)
    
    # Calculate the lambda value for each chamber pressure at each area ratio value
    lambda_1_bar_j = 1 / (m.sqrt(2) * n_1_bar_j * m.pi * d_1_bar_j**2)
    lambda_10_bar_j = 1 / (m.sqrt(2) * n_10_bar_j * m.pi * d_10_bar_j**2)
    lambda_100_bar_j = 1 / (m.sqrt(2) * n_100_bar_j * m.pi * d_100_bar_j**2)
    
    # Append the lambda value to the corresponding lambda array
    lambda_1_bar.append(lambda_1_bar_j)
    lambda_10_bar.append(lambda_10_bar_j)
    lambda_100_bar.append(lambda_100_bar_j)
    
    # Calculate the diameters for each chamber pressure along the nozzle
    nozzle_dimater_1_bar.append(x_to_Diameter(At_vals[0],x_values_1bar[j],theta))
    nozzle_dimater_10_bar.append(x_to_Diameter(At_vals[1],x_values_1bar[j],theta))
    nozzle_dimater_100_bar.append(x_to_Diameter(At_vals[2],x_values_1bar[j],theta))
    
    # Calculate the mean thermal speed of a molecule for each chamber pressure at different
    # area ratios. Units of (m/s)
    mean_therm_speed_1_bar.append(np.sqrt((8*k_b*T_1_bar[j])/(m.pi*mw[0])))
    mean_therm_speed_10_bar.append(np.sqrt((8*k_b*T_10_bar[j])/(m.pi*mw[1])))
    mean_therm_speed_100_bar.append(np.sqrt((8*k_b*T_100_bar[j])/(m.pi*mw[2])))

    # Calculate the local Knudsen number
    knudsen_local_1_bar.append(lambda_1_bar_j / nozzle_dimater_1_bar[j])
    knudsen_local_10_bar.append(lambda_10_bar_j / nozzle_dimater_10_bar[j])
    knudsen_local_100_bar.append(lambda_100_bar_j / nozzle_dimater_100_bar[j])
    
    # Calculate the local collision rate using the formula ... c'(x) / lambda(x)
    collision_rate_local_1_bar.append(mean_therm_speed_1_bar[j] / knudsen_local_1_bar[j])
    collision_rate_local_10_bar.append(mean_therm_speed_10_bar[j] / knudsen_local_10_bar[j])
    collision_rate_local_100_bar.append(mean_therm_speed_100_bar[j] / knudsen_local_100_bar[j])
    
    # Calculate the cummulative number of collisions undergone by a molecule through the nozzle.
    # Use an if-statement to jump the first iteration, since we initialize with zero.
    if j != 0:
        
        # Use the trapezoidal rule to estimate the area under the points
        
        # Calculate the trapezoidal area for chamber pressure 1 bar
        y_val_2 = collision_rate_local_1_bar[j] / u_1_bar[j]
        y_val_1 = collision_rate_local_1_bar[j-1] / u_1_bar[j-1]
        cummulative_collisions_1_bar.append((0.5*del_x_1_bar*(y_val_2+y_val_1)) + cummulative_collisions_1_bar[j-1])
        
        # Calculate the trapezoidal area for chamber pressure 1 bar
        y_val_2 = collision_rate_local_10_bar[j] / u_10_bar[j]
        y_val_1 = collision_rate_local_10_bar[j-1] / u_10_bar[j-1]
        cummulative_collisions_10_bar.append((0.5*del_x_10_bar*(y_val_2+y_val_1)) + cummulative_collisions_10_bar[j-1])
        
        # Calculate the trapezoidal area for chamber pressure 1 bar
        y_val_2 = collision_rate_local_100_bar[j] / u_100_bar[j]
        y_val_1 = collision_rate_local_100_bar[j-1] / u_100_bar[j-1]
        cummulative_collisions_100_bar.append((0.5*del_x_100_bar*(y_val_2+y_val_1)) + cummulative_collisions_100_bar[j-1])



# Define the lengths corresponding to an area ratio of 1.1 from above
x_value_beginner_1_bar   = Area_ratio_to_x(1.1,At_vals[0],theta)
x_value_beginner_10_bar  = Area_ratio_to_x(1.1,At_vals[1],theta)
x_value_beginner_100_bar = Area_ratio_to_x(1.1,At_vals[2],theta)


# Find out at what x values correspond to an A(x)/At value of 1.1
print("\n")
print("An area ratio of 1.1 corresponds to an x value of " + str(x_value_beginner_1_bar) + "m for 1 bar")
print("An area ratio of 1.1 corresponds to an x value of " + str(x_value_beginner_10_bar) + "m for 10 bar")
print("An area ratio of 1.1 corresponds to an x value of " + str(x_value_beginner_100_bar) + "m for 100 bar")


# Find the first element's index of each chamber pressure that is greater than the x values returned
# from the previous step. We have to start from the 9th index for every chamber pressure.
index_beginner_1_bar = np.argmax(x_values_1bar > x_value_beginner_1_bar)
index_beginner_10_bar = np.argmax(x_values_10bar > x_value_beginner_10_bar)
index_beginner_100_bar = np.argmax(x_values_100bar > x_value_beginner_100_bar)
index = 9


# Print out the index
"""
print("Indices are ... ")
print(index_beginner_1_bar)
print(index_beginner_10_bar)
print(index_beginner_100_bar)
"""



# Since the only diatomic species with a molar concentration greater than 1%
# in the product state for all 3 chamber pressures is H2, we will be considering
# H2 collisions with other H2 molecules only.



# Initialize the expansion and relaxation time arrays
exp_time_1_bar   = []
exp_time_10_bar  = []
exp_time_100_bar = []

exp_time_break_down_1_bar   = []
exp_time_break_down_10_bar  = []
exp_time_break_down_100_bar = []

vib_relax_time_1_bar   = []
vib_relax_time_10_bar  = []
vib_relax_time_100_bar = []

rot_relax_time_1_bar   = []
rot_relax_time_10_bar  = []
rot_relax_time_100_bar = []

trans_relax_time_1_bar   = []
trans_relax_time_10_bar  = []
trans_relax_time_100_bar = []


# Initialize important variables to calculate the relaxation time for vibration
charact_vib_temp = 6332                                                             # K
mass_H2          = 1.008                                                            # g / mole
reduced_mass_H2  = mass_H2 / 2                                                      # g / mole
b_AB             = 0.015*((reduced_mass_H2)**(1/4))                                 
a_AB             = 1.16e-3*np.sqrt(reduced_mass_H2)*((charact_vib_temp)**(4/3))



# Pad the arrays with zeros to account for the index offset
for i in range(0,index):
    
    # Append the zeros
    exp_time_1_bar.append(0)
    exp_time_10_bar.append(0)
    exp_time_100_bar.append(0)
    
    exp_time_break_down_1_bar.append(0)
    exp_time_break_down_10_bar.append(0)
    exp_time_break_down_100_bar.append(0)
    
    vib_relax_time_1_bar.append(0)
    vib_relax_time_10_bar.append(0)
    vib_relax_time_100_bar.append(0)
    
    trans_relax_time_1_bar.append(0)
    trans_relax_time_10_bar.append(0)
    trans_relax_time_100_bar.append(0)
    
    rot_relax_time_1_bar.append(0)
    rot_relax_time_10_bar.append(0)
    rot_relax_time_100_bar.append(0)
    
    


# Run a for loop through the diverging nozzle (Starting from index 9)
for i in range(index,len(x_values_1bar)):
    
    # Use an if-statement to determine if we are at the exit area boundary
    if i == (len(x_values_1bar) - 1):
        
        # Use backward finite difference
        exp_time_1_bar.append(exp_time_funct(rho_1_bar[i-1],rho_1_bar[i],del_x_1_bar,u_1_bar[i]))
        exp_time_10_bar.append(exp_time_funct(rho_10_bar[i-1],rho_10_bar[i],del_x_10_bar,u_10_bar[i]))
        exp_time_100_bar.append(exp_time_funct(rho_100_bar[i-1],rho_100_bar[i],del_x_100_bar,u_100_bar[i]))
        
    else:
    
        # Use forward finite difference
        exp_time_1_bar.append(exp_time_funct(rho_1_bar[i],rho_1_bar[i+1],del_x_1_bar,u_1_bar[i]))
        exp_time_10_bar.append(exp_time_funct(rho_10_bar[i],rho_10_bar[i+1],del_x_10_bar,u_10_bar[i]))
        exp_time_100_bar.append(exp_time_funct(rho_100_bar[i],rho_100_bar[i+1],del_x_100_bar,u_100_bar[i]))
    
    
    # Append the corresponding equilibrium breakdown values to the appropriate arrays
    exp_time_break_down_1_bar.append(exp_time_1_bar[i] / 25)
    exp_time_break_down_10_bar.append(exp_time_10_bar[i] / 25)
    exp_time_break_down_100_bar.append(exp_time_100_bar[i] / 25)
    
    # Append the corresponding vibrational relaxation time values to the appropriate arrays 
    vib_relax_time_1_bar.append(relax_time_vib(a_AB,b_AB,P_1_bar[i]*Pa_to_atm,T_1_bar[i]))
    vib_relax_time_10_bar.append(relax_time_vib(a_AB,b_AB,P_10_bar[i]*Pa_to_atm,T_10_bar[i]))
    vib_relax_time_100_bar.append(relax_time_vib(a_AB,b_AB,P_100_bar[i]*Pa_to_atm,T_100_bar[i]))
    
    # Append the corresponding translational relaxation time values to the appropriate arrays 
    trans_relax_time_1_bar.append(1 / collision_rate_local_1_bar[i])
    trans_relax_time_10_bar.append(1 / collision_rate_local_10_bar[i])
    trans_relax_time_100_bar.append(1 / collision_rate_local_100_bar[i])
    
    # Append the corresponding rotational relaxation time values to the appropriate arrays 
    rot_relax_time_1_bar.append(5 / collision_rate_local_1_bar[i])
    rot_relax_time_10_bar.append(5 / collision_rate_local_10_bar[i])
    rot_relax_time_100_bar.append(5 / collision_rate_local_100_bar[i])


################################################ Plotting ###############################################




left_x = [
    [x_vals_nondimen_1_bar, x_vals_nondimen_10_bar, x_vals_nondimen_100_bar],  # plot 1 (placeholders)
    [x_vals_nondimen_1_bar, x_vals_nondimen_10_bar, x_vals_nondimen_100_bar],  # plot 2 (placeholders)
    [x_vals_nondimen_1_bar, x_vals_nondimen_10_bar, x_vals_nondimen_100_bar],  # plot 3 (placeholders)
]
 
# >>> PLACE YOUR LEFT-COLUMN Y ARRAYS HERE (3 plots x 3 curves) <<<
# Each y-array must have the same length as its matching x-array above.
left_y = [
    [knudsen_local_1_bar,knudsen_local_10_bar,knudsen_local_100_bar],              # plot 1 (placeholders)
    [collision_rate_local_1_bar,collision_rate_local_10_bar,collision_rate_local_100_bar],               # plot 2 (placeholders)
    [cummulative_collisions_1_bar,cummulative_collisions_10_bar,cummulative_collisions_100_bar],  # plot 3 (placeholders)
]
left_curve_labels = ["1 bar", "10 bar", "100 bar"]          # shared legend entries
left_curve_colors = ["tab:blue", "tab:orange", "tab:green"]
 
left_title  = "Diverging Nozzle 50 mN Design"
left_xlabel = r"Axial distance from throat $x/D_t$"
left_ylabels = ["Local Kn", r"Collision rate per molecule $\nu$ (s$^{-1}$)",
                "Cumulative collisions\nper molecule from throat"]
 
# Optional fixed y-limits per plot: (ymin, ymax) or None for autoscale
left_ylims = [None,None,None]
# X-limits shared by the 3 left plots: (xmin, xmax), or None = min/max over ALL left x data
left_xlim = None
 
# Optional shaded horizontal bands per plot: list of (ymin, ymax, color, alpha)
left_bands = [
    [(1e-1, 1e0, "gray", 0.69), (1e-2, 1e-1, "gray", 0.15)],   # plot 1 (example)
    [],                                                        # plot 2
    [],                                                        # plot 3
]
 
# ---------------------------------------------------------------------
# RIGHT COLUMN DATA
# ---------------------------------------------------------------------
# Each curve has its OWN x-array and y-array. The plots only share the x-axis
# LIMITS/ticks (see right_xlim below), not the x data itself.
# right_x[plot_index][curve_index] pairs with right_y[plot_index][curve_index]
# Curve order must match right_curve_styles below. Y values must be positive.
# >>> PLACE YOUR RIGHT-COLUMN X ARRAYS HERE (3 plots x 5 curves) <<<
def _x(start, stop, n):
    return np.linspace(start, stop, n)
right_x = [
    [x_vals_nondimen_1_bar,x_vals_nondimen_1_bar,x_vals_nondimen_1_bar,x_vals_nondimen_1_bar,x_vals_nondimen_1_bar],  # plot 1 (placeholders)
    [x_vals_nondimen_10_bar,x_vals_nondimen_10_bar,x_vals_nondimen_10_bar,x_vals_nondimen_10_bar,x_vals_nondimen_10_bar],  # plot 2 (placeholders)
    [x_vals_nondimen_100_bar,x_vals_nondimen_100_bar,x_vals_nondimen_100_bar,x_vals_nondimen_100_bar,x_vals_nondimen_100_bar],  # plot 3 (placeholders)
]
 
# >>> PLACE YOUR RIGHT-COLUMN Y ARRAYS HERE (3 plots x 5 curves) <<<
# Each y-array must have the same length as its matching x-array above.
_scales = [1.0, 0.5, 0.1]
_factors = [1, 1/25, 1/80, 1/10, 1e2]      # placeholder curve shapes
right_y = [
    [exp_time_1_bar,exp_time_break_down_1_bar,vib_relax_time_1_bar,trans_relax_time_1_bar,rot_relax_time_1_bar],  # plot 1 (placeholders)
    [exp_time_10_bar,exp_time_break_down_10_bar,vib_relax_time_10_bar,trans_relax_time_10_bar,rot_relax_time_10_bar],  # plot 2 (placeholders)
    [exp_time_100_bar,exp_time_break_down_100_bar,vib_relax_time_100_bar,trans_relax_time_100_bar,rot_relax_time_100_bar],  # plot 3 (placeholders)
]


 
# One entry per curve (label, color, linestyle, linewidth)
right_curve_styles = [
    dict(label=r"$\tau_{exp}$ (time available)", color="k",          ls="-",  lw=2.5),
    dict(label=r"$\tau_{exp}/25$",               color="k",          ls="--", lw=1.2),
    dict(label=r"$\tau_{v}\approx 1/\nu$",   color="dimgray",    ls="-",  lw=1.5),
    dict(label=r"$\tau_{trans}\approx 5/\nu$",     color="tab:orange", ls="-",  lw=1.5),
    dict(label=r"$\tau_{rot}$ (CO)",               color="tab:red",    ls="--", lw=1.5),
]
 
# Optional shaded regions between curves (same for every right plot).
# (lower, upper, color, alpha, label); lower/upper = curve index,
# or "bottom" / "top" for the edge of the plot. Use [] for none.
right_regions = [
    ("bottom", 1,     "tab:green",  0.15, r"Equilibrium ($\tau_k<\tau_{exp}/25$)"),
    (1,        0,     "gold",       0.30, r"Breakdown ($\tau_{exp}/25<\tau_k<\tau_{exp}$)"),
    (0,        "top", "slategray",  0.20, r"Frozen ($\tau_k>\tau_{exp}$)"),
]
 
right_title  = "Time Needed to Relax vs Time Available 50 mN Design"
right_xlabel = r"Axial distance from throat $x/D_t$"
right_ylabels = ["Time (s)", "Time (s)", "Time (s)"]
 
# >>> TYPE YOUR TEXT BOX CONTENTS HERE (one string per right plot) <<<
right_textbox_text = ["1 bar", "10 bar", "100 bar"]
 
# >>> SET THE Y-LIMITS SHARED BY ALL THREE RIGHT PLOTS HERE <<<
# (ymin, ymax), both positive (log axis). Same range on every right plot so
# they are directly comparable. Set to None to autoscale each plot instead.
right_ylim = (1e-10, 1e-3)
# X-limits shared by the 3 right plots: (xmin, xmax), or None = min/max over ALL right x data
right_xlim = None
 
 
# =====================================================================
# HELPERS
# =====================================================================
def style_axes(ax, ylabel, ylim):
    ax.set_yscale("log")
    ax.grid(True, which="major", alpha=0.35)
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(ylim)
 
 
def log_interp(xg, x, y):
    """Interpolate y(x) onto xg in log space (y must be positive)."""
    return 10**np.interp(xg, x, np.log10(y))
 
 
def fill_between_curves(ax, xs, ys, lo_idx, hi_idx, lo, hi, **kw):
    """Shade between two curves that may have different x-arrays.
    lo_idx / hi_idx: curve index, or "bottom" / "top" for the plot edge."""
    idx = [k for k in (lo_idx, hi_idx) if k not in ("bottom", "top")]
    start = max(xs[k][0] for k in idx)       # overlap of the x-ranges involved
    stop = min(xs[k][-1] for k in idx)
    xg = np.unique(np.concatenate([xs[k][(xs[k] >= start) & (xs[k] <= stop)] for k in idx]))
    y_lo = np.full_like(xg, lo) if lo_idx == "bottom" else log_interp(xg, xs[lo_idx], ys[lo_idx])
    y_hi = np.full_like(xg, hi) if hi_idx == "top" else log_interp(xg, xs[hi_idx], ys[hi_idx])
    ax.fill_between(xg, y_lo, y_hi, **kw)
 
 
def add_textbox(ax, text):
    ax.text(0.98, 0.05, text, transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="square", facecolor="white", edgecolor="0.7"))
 
 
def finish_column(axes, x_arrays, xlim, xlabel, title):
    """x_arrays: every x-array used in the column (used only for the default limits)."""
    if xlim is None:
        allx = np.concatenate([np.ravel(a) for a in x_arrays])
        xlim = (np.min(allx), np.max(allx))
    axes[0].set_xlim(xlim)                       # sharex applies it to all three
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)        # tick marks stay, numbers hidden
    axes[-1].set_xlabel(xlabel)
    axes[0].set_title(title, fontsize=13)
 
 
# =====================================================================
# BUILD THE FIGURE
# =====================================================================
fig, axs = plt.subplots(3, 2, sharex="col", figsize=(16, 11))
axL, axR = axs[:, 0], axs[:, 1]
 
# ----- LEFT column -----------------------------------------------------
for i, ax in enumerate(axL):
    for x, y, lab, col in zip(left_x[i], left_y[i], left_curve_labels, left_curve_colors):
        ax.plot(x, y, label=lab, color=col, zorder=3)
    style_axes(ax, left_ylabels[i], left_ylims[i])
    lo, hi = ax.get_ylim()
    for (y0, y1, c, a) in left_bands[i]:
        ax.axhspan(y0, y1, color=c, alpha=a, lw=0, zorder=0)
    ax.set_ylim(lo, hi)                          # bands must not change the limits
finish_column(axL, [x for row in left_x for x in row], left_xlim, left_xlabel, left_title)
axL[0].legend(loc="upper left", fontsize=8)      # shared legend, inside top plot
 
# ----- RIGHT column ----------------------------------------------------
for i, ax in enumerate(axR):
    for x, y, st in zip(right_x[i], right_y[i], right_curve_styles):
        ax.plot(x, y, color=st["color"], ls=st["ls"], lw=st["lw"],
                label=st["label"], zorder=3)
    style_axes(ax, right_ylabels[i], right_ylim)
    lo, hi = ax.get_ylim()
    for (a, b, c, alpha, lab) in right_regions:
        fill_between_curves(ax, right_x[i], right_y[i], a, b, lo, hi,
                            color=c, alpha=alpha, lw=0,
                            label=lab if i == 0 else None, zorder=0)
    ax.set_ylim(lo, hi)
    add_textbox(ax, right_textbox_text[i])
finish_column(axR, [x for row in right_x for x in row], right_xlim, right_xlabel, right_title)
 
# Shared legend for the right column, centered under the bottom plot
handles, labels = axR[0].get_legend_handles_labels()
axR[-1].legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.28),
               ncol=2, fontsize=8)
 
fig.subplots_adjust(left=0.06, right=0.98, top=0.94, bottom=0.20,
                    wspace=0.25, hspace=0.12)
plt.show()
