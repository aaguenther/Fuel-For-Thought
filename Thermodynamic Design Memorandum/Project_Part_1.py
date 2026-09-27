import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

######################################################FUNCTIONS#########################################################

# Define the constant pressure adiabatic flame temperature function
def equilibrium_properties(fuel, oxidizer, T, P, eq_ratio, gas_object):
    
    # Initialize arrays for AFT, gamma, and molecular weight
    flame_temp = []
    gamma = []
    mw = []

    # Set the initial temperature and pressure for the reactants
    gas_object.TP = T, P

    # Run a for loop for every equivalence ratio value and compute the equilibrium conditions
    for i in range(0,len(eq_ratio)):

        # Reset the object state after each iteration
        gas_object.TP = T, P

        # Set up the equivalence ratio
        gas_object.set_equivalence_ratio(phi=eq_ratio[i], fuel=fuel, oxidizer=oxidizer)
        
        # Find the equilibrium of the reaction to find the flame temperature
        gas_object.equilibrate('HP')

        # Store the calculated properties in arrays
        flame_temp.append(gas_object.T)
        gamma.append(gas_object.cp_mass / gas_object.cv_mass)
        mw.append(gas_object.mean_molecular_weight) #kg/kmol

    return flame_temp, gamma, mw

# Define the code to create the plot for the adiabatic flame temperatures
def plot_flame_temp_curves(eq_ratio,flame_temp_curves,P,title):

    # Set up the plotting figure
    plt.figure(figsize=(10,8))

    # Plot the three different flame temp curves for each pressure
    plt.plot(eq_ratio,flame_temp_curves[0][0],label=str(P[0]/1E5)+' Bar',color='blue',linestyle='-')
    plt.plot(eq_ratio,flame_temp_curves[1][0],label=str(P[1]/1E5)+' Bar',color='purple',linestyle='-')
    plt.plot(eq_ratio,flame_temp_curves[2][0],label=str(P[2]/1E5)+' Bar',color='red',linestyle='-')

    # Format the layout of the plot
    plt.title(title)
    plt.xlabel('Equivalence Ratio (phi)')
    plt.ylabel('Adiabatic Flame Temperature (K)')
    plt.xlim(min(eq_ratio)-0.2, max(eq_ratio)+0.2)
    plt.grid(True)
    plt.legend()

    plt.show()

# Function to calculate C*
def characteristic_velocity(properties):

    # From the tuple of properties, create arrays for each property
    T_ad = np.array(properties[0])
    gamma = np.array(properties[1])
    mw = np.array(properties[2])
    R = Ru/mw

    # Calculate C*
    c_star = (1 / gamma) * ((gamma + 1) / 2) ** ((gamma + 1) / (2 * (gamma - 1))) * (gamma * R * T_ad) ** 0.5

    return c_star

# Function to plot C*
def plot_characteristic_velocity(eq_ratio, c_star, P, title):

    # Initialize figure
    plt.figure(figsize=(10,8))

    # Plot c_star for each pressure

    plt.plot(eq_ratio, c_star[0], label=str(P[0] / 1E5) + ' Bar', color='blue', linestyle='-')
    plt.plot(eq_ratio, c_star[1], label=str(P[1] / 1E5) + ' Bar', color='purple', linestyle='-')
    plt.plot(eq_ratio, c_star[2], label=str(P[2] / 1E5) + ' Bar', color='red', linestyle='-')

    # Format the layout of the plot
    plt.title(title)
    plt.xlabel('Equivalence Ratio (phi)')
    plt.ylabel('c* (m/s)')
    plt.xlim(min(eq_ratio) - 0.2, max(eq_ratio) + 0.2)
    plt.grid(True)
    plt.legend()

    plt.show()

#################################################CALCULATE AND PLOT#####################################################

# Create two objects
gas_1 = ct.Solution("gri30.yaml") # Mechanism for H2 and CH4
gas_2 = ct.Solution("nDodecane_Reitz.yaml") # Mechanism for kerosene

# Populate the equivalence ratio array
eq_ratio = np.linspace(0.5,4, 60)

# Chamber temperature and pressure
pressures = [0.5E7,1.5E7,3E7] # (Pa)
temperature = 298 # (K)
Ru = 8314 # J/(kmol*K)

# Define the different fuel and oxidizer strings
kerosene_2 = "c12h26"
methane_1  = "CH4"
hydrogen_1 = "H2"
oxygen_2   = "o2" # Oxygen is named differently in the two mechanisms
oxygen_1   = "O2"

# Initialize the array of properties for kerosene, methane, and hydrogen
kerosene_properties = []
methane_properties = []
hydrogen_properties = []

# Populate the arrays
for i in range(0,len(pressures)):
    
    # Append the equilibrium properties to the arrays
    kerosene_properties.append(equilibrium_properties(kerosene_2,oxygen_2,temperature,pressures[i],eq_ratio,gas_2))
    methane_properties.append(equilibrium_properties(methane_1,oxygen_1,temperature,pressures[i],eq_ratio,gas_1))
    hydrogen_properties.append(equilibrium_properties(hydrogen_1,oxygen_1,temperature,pressures[i],eq_ratio,gas_1))
    # This creates a list of tuples, where each tuple contains three arrays corresponding to one pressure

# Define arrays of characteristic velocities for different fuels
kerosene_c_star = []
methane_c_star = []
hydrogen_c_star = []

# Populate the arrays
for i in range(0, len(pressures)):
    # Append c_star to arrays
    kerosene_c_star.append(characteristic_velocity(kerosene_properties[i]))
    methane_c_star.append(characteristic_velocity(methane_properties[i]))
    hydrogen_c_star.append(characteristic_velocity(hydrogen_properties[i]))
    # This creates a list of arrays, where each corresponds to one pressure

# Define the title strings for AFT
kerosene_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Kerosene + Oxygen)"
methane_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Methane + Oxygen)"
hydrogen_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Hydrogen + Oxygen)"

# Define the title strings for c*
kerosene_c_star_title = "Characteristic Velocity vs. Equivalence Ratio (Kerosene + Oxygen)"
methane_c_star_title = "Characteristic Velocity vs. Equivalence Ratio (Methane + Oxygen)"
hydrogen_c_star_title = "Characteristic Velocity vs. Equivalence Ratio (Hydrogen + Oxygen)"

# Call the plotting function to plot flame temp curves versus equivalence ratio
plot_flame_temp_curves(eq_ratio,kerosene_properties,pressures,kerosene_title)
plot_flame_temp_curves(eq_ratio,methane_properties,pressures,methane_title)
plot_flame_temp_curves(eq_ratio,hydrogen_properties,pressures,hydrogen_title)

# Call plotting function to plot c* vs equivalence ratio
plot_characteristic_velocity(eq_ratio,kerosene_c_star, pressures, kerosene_c_star_title)
plot_characteristic_velocity(eq_ratio,methane_c_star, pressures, methane_c_star_title)
plot_characteristic_velocity(eq_ratio,hydrogen_c_star, pressures, hydrogen_c_star_title)



