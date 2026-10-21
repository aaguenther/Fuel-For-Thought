import math as m
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Given

T_initial = 298 # K
P_c_bar = np.array([1, 10, 100]) # bar
P_c = P_c_bar * 1e5 # Pa
phi = [3.0]
R_u = ct.gas_constant # J / (kmol * K)
N_A = 6.0221408e23
k_b = 1.38064852e-23

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