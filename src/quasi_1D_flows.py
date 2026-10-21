from scipy.optimize import brentq

def area_mach_relation(area_ratio, gamma, branch):
    """Function to compute either supersonic or subsonic mach number solution for a given area ratio and gamma"""

    def f(mach):
        return area_ratio - (1 / mach) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * mach**2)) ** ((gamma + 1) / (2 * (gamma - 1)))
    # Sonic condition
    if abs(area_ratio - 1.0) < 1e-10:
        return 1.0
    if branch == "subsonic":
        M = brentq(f, 1E-6, 0.999)
    elif branch == "supersonic":
        M = brentq(f, 1.001, 10)
    else:
        print("branch must be 'supersonic' or 'subsonic'")
        return None

    return M


