def po_p (gamma, m):
    """Calculate pressure ratio Po/P for an isentropic flow"""
    return (1 + 0.5 * (gamma - 1) * m**2)**(gamma / (gamma - 1))

def speed_of_sound (gamma, r, t):
    """Calculate Mach number for an isentropic flow"""
    return (gamma * r * t)**0.5

def stagnation_temperature (t, u, cp):
    """Calculate stagnation temperature for an isentropic flow"""
    return t + u**2 / (2 * cp)

def to_t (gamma, m):
    """Calculate temperature ratio To/T for an isentropic flow"""
    return 1 + 0.5 * (gamma - 1) * m**2

def rho_o_rho (gamma, m):
    """Calculate density ratio rho_o/rho for an isentropic flow"""
    return (1 + 0.5 * (gamma - 1) * m**2) ** (1 / (gamma - 1))

def t2_t1 (gamma, m1, m2):
    """Calculate temperature ratio T2/T1 for an isentropic flow"""
    return (1 + 0.5 * (gamma - 1) * m1**2) / (1 + 0.5 * (gamma - 1) * m2**2)

def p2_p1 (gamma, m1, m2):
    """Calculate pressure ratio P2/P1 for an isentropic flow"""
    return t2_t1(gamma, m1, m2) ** (gamma / (gamma - 1))

def F (gamma, M):
    """Calculate area ratio for an isentropic flow"""
    return (1 / M) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M**2)) ** ((gamma + 1) / (2 * (gamma - 1)))

