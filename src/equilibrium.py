def equilibrium_properties(fuel, oxidizer, T, P, eq_ratio, gas_object):
    # Initialize arrays for AFT, gamma, and molecular weight
    flame_temp = []
    gamma = []
    mw = []

    # Set the initial temperature and pressure for the reactants
    gas_object.TP = T, P

    # Run a for loop for every equivalence ratio value and compute the equilibrium conditions
    for i in range(0, len(eq_ratio)):
        # Reset the object state after each iteration
        gas_object.TP = T, P

        # Set up the equivalence ratio
        gas_object.set_equivalence_ratio(phi=eq_ratio[i], fuel=fuel, oxidizer=oxidizer)

        # Find the equilibrium of the reaction to find the flame temperature
        gas_object.equilibrate('HP')

        # Store the calculated properties in arrays
        flame_temp.append(gas_object.T)
        gamma.append(gas_object.cp_mass / gas_object.cv_mass)
        mw.append(gas_object.mean_molecular_weight)  # kg/kmol

    return flame_temp, gamma, mw, gas_object

def characteristic_velocity(T_ad, gamma, R):

    # Calculate C*
    c_star = (1 / gamma) * ((gamma + 1) / 2) ** ((gamma + 1) / (2 * (gamma - 1))) * (gamma * R * T_ad) ** 0.5

    return c_star