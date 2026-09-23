import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


gas_1 = ct.Solution("gri30.yaml")
gas_2 = ct.Solution("nDodecane_Reitz.yaml")


"""
print(gas_1.species_names)
print("\n\n\n")
print(gas_2.species_names)
"""


# The main goal of this script is to calculate the constant pressure adiabatic flame temperatures
# for three different combustion reactions (kerosine, methane, and hydrogen).


# The second goal of this script is to create 6 different plots: 3 regarding the flame temperatures
# of every reaction as a function of the equivalence ratio and 3 regarding the c_star speed
# of every reaction as a function of the equivalence ratio. Each graph will have three different curves
# corresponding to different pressure levels


# ________________________________________________________________________________________________________

# Create a function called constant_pressure_flame_temperature() with the following inputs and outputs
# ________________________________________________________________________________________________________
#
#                    Inputs                                                 Outputs
#           ________________________                 |             __________________________
#                                                    |
#   A list containing strings of the chemical        |          A list of the flame temperatures as 
#   species and their stoich. coefficients for       |          a function of the different equivalence
#   a given reaction.                                |          ratio values
#                                                    |
#   A list of the different equivalence ratios       |
#   that will be used for analysis                   |
#                                                    |
#   A list containing the numeric values for the     |
#   initial temperature and pressure of the          |
#   reaction                                         |
# ________________________________________________________________________________________________________

# Set up important parameters and variables for the analyses and plots
# ________________________________________________________________________________________________________

# Populate the equivalence ratios
eq_ratio = np.linspace(0.3,3,5)


# Populate the different chamber pressures to use in [Pa] and the chamber temperature [K]
pressures = [101325,500000,1000000]
temperature = 298


# Populate the different fuel and oxidizer strings
kerosine_2 = "c12h26"
methane_1  = "CH4"
hydrogen_1 = "H2"
oxygen_2   = "o2"
oxygen_1   = "O2"

# ________________________________________________________________________________________________________

# Define important functions to help carry out the analyses and visual plotting
# ________________________________________________________________________________________________________


# Define the constant pressure adiabatic flame temperature function
def constant_pressure_flame_temperature(fuel,oxidizer,T,P,eq_ratio,gas_object):
    
    
    # Initialize the adiabatic flame temperature array that will be returned
    flame_temp = []
    
    
    # Set the initial temperature and pressure for the reactants
    gas_object.TP = T, P
    
    
    # Run a for loop for every equivalence ratio value and compute the adiabatic flame temperature
    for i in range(0,len(eq_ratio)):
        
        
        # Set up the equivalence ratio
        gas_object.set_equivalence_ratio(phi=eq_ratio[i], fuel=fuel, oxidizer=oxidizer)
        
        
        # Find the equilibrium of the reaction to find the flame temperature
        gas_object.equilibrate('HP')
        
        
        # Store the calculated temperature in the temperature array
        flame_temp.append(gas_object.T)
        
        
    return flame_temp


# Define the code to create the plot for the adiabatic flame temperatures
def plot_flame_temp_curves(eq_ratio,flame_temp_curves,P,title):
    
    
    # Set up the plotting figure
    plt.figure(figsize=(10,8))
    
    
    # Plot the three different flame temp curves for each pressure
    plt.plot(eq_ratio,flame_temp_curves[0],label=str(P[0])+' Pa',color='blue',linestyle='-')
    plt.plot(eq_ratio,flame_temp_curves[1],label=str(P[1])+' Pa',color='purple',linestyle='-')
    plt.plot(eq_ratio,flame_temp_curves[2],label=str(P[2])+' Pa',color='red',linestyle='-')
    
    
    # Format the layout of the plot
    plt.title(title)
    plt.xlabel('Equivalence Ratio (phi)')
    plt.ylabel('Flame Temperature (K)')
    plt.xlim(min(eq_ratio), max(eq_ratio))
    plt.grid(True)
    plt.legend()
    
    # Save and show
    #plt.savefig(filename)
    plt.show()



# ________________________________________________________________________________________________________

# Define the different chemical reactions and call the user-defined functions above 
# ________________________________________________________________________________________________________



# Define the adiabatic flame temperature array for kerosine, methane, and hydrogen
kerosine_flame_temp_array = []
methane_flame_temp_array = []
hydrogen_flame_temp_array = []

# Populate the array
for i in range(0,len(pressures)):
    
    # Append the flame temp outputs into the kerosine, methane, and hydrogen array
    kerosine_flame_temp_array.append(constant_pressure_flame_temperature(kerosine_2,oxygen_2,temperature,pressures[i],eq_ratio,gas_2))
    methane_flame_temp_array.append(constant_pressure_flame_temperature(methane_1,oxygen_1,temperature,pressures[i],eq_ratio,gas_1))
    hydrogen_flame_temp_array.append(constant_pressure_flame_temperature(hydrogen_1,oxygen_1,temperature,pressures[i],eq_ratio,gas_1))


# Define the title strings
kerosine_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Kerosine + Oxygen)"
methane_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Methane + Oxygen)"
hydrogen_title = "Adiabatic Flame Temperature vs. Equivalence Ratio (Hydrogen + Oxygen)"


# Call the plotting function to print out the flame temp curves versus equivalence ratio
plot_flame_temp_curves(eq_ratio,kerosine_flame_temp_array,pressures,kerosine_title)
plot_flame_temp_curves(eq_ratio,methane_flame_temp_array,pressures,methane_title)
plot_flame_temp_curves(eq_ratio,hydrogen_flame_temp_array,pressures,hydrogen_title)