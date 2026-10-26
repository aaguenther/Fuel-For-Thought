import numpy as np
import math as m


def x_to_Area_Ratio(At,x,theta):
    
    # This function will take in an array of x values (length along the nozzle length) and return
    # an array containing the corresponding cross-sectional area ratio (A(x)/At) values at that location
    
    # The inputs are the throat area, the x values array, and the angle theta of the diverging nozzle.
    #The x values array must be a numpy array. Theta must be in radians.
    
    # Apply the formula for the corresponding formula for the area ratio
    
    # Obtain Dt
    Dt = np.sqrt(4*At/m.pi)
    
    # Calculate the first term
    A_ratio_1 = ((2*np.tan(theta)/Dt)*x)**2
    
    # Calculate the second term
    A_ratio_2 = (4*np.tan(theta)/Dt)*x
    
    # Third term is a constant
    A_ratio_3 = 1
    
    # Bring all of the terms together
    A_ratio = A_ratio_1 + A_ratio_2 + A_ratio_3
    
    # Return the area ratios array
    return A_ratio




def nozzle_length(Ae,At,theta):
    
    # This function will return a scale value representing the length of the diverging nozzle
    # section
    
    # The inputs are both scalar: the exit and throat area.
    
    # Apply the formula
    L = (np.sqrt(4*Ae/m.pi) - np.sqrt(4*At/m.pi)) / (2*np.tan(theta))
    
    # Returnt the result
    return L




def x_to_Diameter(At,x,theta):
    
    # This function will convert x values into their corresponding diameter values
    # The inputs needed are the throat area, x value along the nozzle, and the theta value
    
    # Apply the formula for the diameter
    D = np.sqrt(4*At/m.pi) + 2*x*np.tan(theta)
    
    # Return the diameter
    return D




def Area_ratio_to_x(ratio,At,theta):
    
    # This function will serve as the inverse of the x_to_Area_Ratio function
    # The quadratic function will be used to find the corresponding x value
    # The ratio input is defined as A(x) / At
    
    # Define the throat diameter
    Dt = np.sqrt(4*At/m.pi)
    
    # Define a
    a = (2*np.tan(theta)/Dt)**2
    
    # Define b
    b = (4*np.tan(theta)/Dt)
    
    # Define c
    c = 1 - ratio
    
    # Define the determinant
    determinant = np.sqrt(b**2 - (4*a*c))
    
    # Implement the quadratic formula to obtain x
    x = (-b + determinant) / (2*a)
    
    # Check to see if x is < 0
    if x < 0:
        print("Something is not right!!!\n")
    
    # Return the x value
    return x





def exp_time_funct(den1,den2,del_x,u):
    
    # This function will return the expansion time value. Can use either forward or backward
    # Finite difference to calculate the derivative term.
    
    # Den2 refers to f(x_i+1) & f(x_i) for forward and backward schemes respectively
    # Den1 refers to f(x_i) & f(x_i-1) for forward and backward schemes respectively
        
    # Implement the forward finite difference simplification
    ans = 1 / np.abs( (u * np.log(den2 / den1)) / (del_x) )
    
    # Return the results
    return ans





def relax_time_vib(a,b,P,T):
    
    # This function will return the estimate of the relaxation time using the Millikan–White correlation
    # The pressure input needs to be in atm, T in K
    
    # Implement the formula
    t = np.exp( a*( ((T)**(-1/3)) - b ) - 18.42 ) / P
    
    # Return the relaxation time
    return t