# Authors: Argo Papadopoulou, Thomas Tsekas, Konstantinos Papalamprou
# Date: 25/5/2026
#
# Introduction to ECE - team 09 - BJT
# Find R_2 using Bisection Method

import numpy as np

######################################################################### Constant Data #########################################################################
# Settings
step = 1
small_step = 0.01
accuracy = 0.01
FIG_PATH = "../assets/V(R_2).png"

# BJT
I_s = 10**-16 # A
beta = 100
V_t = 0.026 # V
V_a = 100 # V
V_all = [[], [], []] # All the solutions of the (V_1, V_2, V_3) system that the Newton method will find

# Source
V_cc = 15 # V

# Resistors
R_L = 2000 # Ohm
R_s = 0.1 * R_L # Ohm
R_1 = 10 * R_L # Ohm

R_2_min = 100  # Ohm
R_2_max = 4000 # Ohm
R_2_all = [] # All the value that the resistor will get

# (iii)
target_V_3 = 7.5000

# System of functions
def F(V, R_2):
    V_1, V_2, V_3 = V

    f1 = (V_1/R_2) + ((V_1 - V_cc)/R_1) + (I_s * (np.exp((V_1 - V_2)/V_t) - 1))
    f2 = (V_2/R_s) - (I_s * (np.exp((V_1 - V_2)/V_t) - 1)) - (beta * I_s * (np.exp((V_1 - V_2)/V_t) - 1) * (1 + ((V_3 - V_2) / V_a)))
    f3 = ((V_3 - V_cc) / R_L) + (beta * I_s * (np.exp((V_1 - V_2)/V_t) - 1) * (1 + ((V_3 - V_2) / V_a)))

    return [f1, f2, f3]

# We can generate the jacobian using the [project]/utils/find_j.py
def J(V, R_2):
    V_1, V_2, V_3 = V

    df1_dV_1 = I_s*np.exp((V_1 - V_2)/V_t)/V_t + 1/R_2 + 1/R_1
    df2_dV_1 = -I_s*np.exp((V_1 - V_2)/V_t)/V_t
    df3_dV_1 = 0

    df1_dV_2 = -I_s*beta*(1 + (-V_2 + V_3)/V_a)*np.exp((V_1 - V_2)/V_t)/V_t - I_s*np.exp((V_1 - V_2)/V_t)/V_t
    df2_dV_2 = I_s*beta*(1 + (-V_2 + V_3)/V_a)*np.exp((V_1 - V_2)/V_t)/V_t + I_s*np.exp((V_1 - V_2)/V_t)/V_t + I_s*beta*(np.exp((V_1 - V_2)/V_t) - 1)/V_a + 1/R_s
    df3_dV_2 = -I_s*beta*(np.exp((V_1 - V_2)/V_t) - 1)/V_a

    df1_dV_3 = I_s*beta*(1 + (-V_2 + V_3)/V_a)*np.exp((V_1 - V_2)/V_t)/V_t
    df2_dV_3 = -I_s*beta*(1 + (-V_2 + V_3)/V_a)*np.exp((V_1 - V_2)/V_t)/V_t - I_s*beta*(np.exp((V_1 - V_2)/V_t) - 1)/V_a
    df3_dV_3 = I_s*beta*(np.exp((V_1 - V_2)/V_t) - 1)/V_a + 1/R_L

    return [[df1_dV_1, df2_dV_1, df3_dV_1],
            [df1_dV_2, df2_dV_2, df3_dV_2],
            [df1_dV_3, df2_dV_3, df3_dV_3]]

############################################################### Implementation of Newtons method (i) #############################################################

# Implementation of Newton method: X_n+1 = X_n - J^-1 * F(X_n) <=> X_n+1 - X_n = - J^-1 * F(X_n) <=> delta = - J^-1 * F(X_n) <=> J * delta = -F(X_n) (1)
def newton_solve(F, x0, args=(), tol=1e-10, max_iter=50):
    x = np.array(x0, dtype=float) # Initial coordinates

    for i in range(max_iter):
        Fx = np.array(F(x, *args), dtype=float)

        # Check if the function converged near zero to stop the process
        if (np.linalg.norm(Fx) < tol):
            # print(f"DEBUG: Converged near zero after {i} iterations")
            return x
        
        Jx = np.array(J(x, *args), dtype=float)
        delta = np.linalg.solve(Jx, -Fx) # Solve the system (1) using LU algorithm

        x += delta

    raise RuntimeError(f"Could not corverge for x0={x0} with settings: max_iter={max_iter} and tol={tol}")

######################################################################### Bisection Method #######################################################################

def bisection(target_V_3, low, high, tol=1e-6):
    (V_1_0, V_2_0, V_3_0) = (0, 0, 0) # (V, V, V) NOTE: Initial V & variable to save the values of V in the i-1 iteration

    for _ in range(100):
        mid = (low + high) / 2
        (V_1_mid, V_2_mid, V_3_mid) = newton_solve(F, x0=[V_1_0, V_2_0, V_3_0], args=(mid, ))
        (V_1_0, V_2_0, V_3_0) = (V_1_mid, V_2_mid, V_3_mid)

        if V_3_mid > target_V_3:
            low = mid

        else:
            high = mid

        if abs(V_3_mid - target_V_3) < tol:
            return mid, V_3_mid
        
    return mid, V_3_mid

################################################################### Apply the Bisection Method ####################################################################

R_2, V_3_mid = bisection(target_V_3, R_2_min, R_2_max)
print(f"V_3 = {V_3_mid:.2f} V => R_2 ~= {R_2:.2f} Ω")