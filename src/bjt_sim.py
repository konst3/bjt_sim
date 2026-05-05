# Authors: Konstantinos Papalamprou
# Date: 25/4/2026
#
# Introduction to ECE - team 09 - BJT
# Main BJT simulation

import numpy as np
# from scipy.optimize import fsolve # Solve a system using Newton method (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html)
import matplotlib.pyplot as plt

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
(V_1_0, V_2_0, V_3_0) = (0, 0, 0) # (V, V, V) NOTE: Initial V & variable to save the values of V in the i-1 iteration
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

###################################################################### Solve the System (i) ######################################################################

print(" R_2 (Ω) | V_1 (V) | V_2 (V) | V_3(V) ")
print("-----------------------------------------")

R_2 = R_2_min
stp = step
R_2_intersect = 0

while (R_2_min <= R_2 <= R_2_max):
    (V_1, V_2, V_3) = newton_solve(F, x0=[V_1_0, V_2_0, V_3_0], args=(R_2, )) # NOTE: Use the fsolve function from SciPy to solve the system using Newton method 
    (V_1_0, V_2_0, V_3_0) = (V_1, V_2, V_3)

    print(f" {R_2:7.2f} | {V_1:7.4f} | {V_2:7.4f} | {V_3:7.4f}", end="")

    # (iii) --> Decrease step near target V_3
    if ((target_V_3 - accuracy) <= V_3 <= (target_V_3 + accuracy)):
        print(" --> decreased step", end="")
        stp = small_step

    else:
        R_2 = int(R_2)
        stp = step

    # NOTE: To avoid floating point error calculate the difference and DO NOT USE ==
    if (abs(V_3 - target_V_3) < 5e-5): 
        print(" --> (iii)", end="")
        R_2_intersect = R_2

    # NOTE: Do we need to change equations for the transistor?
    if (V_3 < V_1): print(" --> Invalid: Transistor Saturation!")
    else: print("")

    # Add data to be plotted
    V_all[0].append(V_1)
    V_all[1].append(V_2)
    V_all[2].append(V_3)
    R_2_all.append(R_2)

    R_2 += stp # Solve the system for next resistance

######################################################################### Plot data (ii) #########################################################################

fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(R_2_all, V_all[0], label="$ V_1 $", color="blue", alpha=0.3)
ax.plot(R_2_all, V_all[1], label="$ V_2 $", color="blue", alpha=0.6)
ax.plot(R_2_all, V_all[2], label="$ V_3 $", color="red")

ax.set_ylim(0, 20)
ax.set_xlim(R_2_min, R_2_max)

ax.set(xlabel = "$ R_2 (Ω) $", ylabel = "$ Voltages (V) $",
       title = "$ V(R_2) $")

# (iii)
ax.axhline(y=target_V_3, color="black", linestyle="--", label=f"$ V_3={target_V_3}V $")

# # Find R_2 from the intersection of x=V_3 and the V_3(R_2) --> Now we decrease the step close to the target V_3
# # NOTE: Find the value of V from the V_3 array that is closest the the target voltage
# #       np.diff()     --> get the difference 
# #       np.argwhere() --> get where the sign changes
# #       .flatten()    --> from 2D to 1D
# diff = np.array(V_all[2]) - target_V_3
# idx = np.argwhere(np.diff(np.sign(diff))).flatten()

# #       np.interp()   --> interpolate from the to values with differences with alternating signs (ex. -0.8, +0.2)
# #                         and their corresponding resistance to find the exact resistance where the difference is 0
# R_2_intersect = np.interp(0, [diff[idx[0]], diff[idx[0]+1]], [R_2_all[idx[0]], R_2_all[idx[0]+1]])

print(f"V_3 = {target_V_3} V => R_2 ~= {R_2_intersect:.2f} Ω")

# Mark the point on the plot
ax.axvline(x=R_2_intersect, color="black", linestyle="--")

ax.plot(R_2_intersect, target_V_3, 'ro', markersize=5, label=f"R_2 = {R_2_intersect:.1f} Ω")

ax.legend(loc="best")
ax.grid()

fig.savefig(FIG_PATH)
plt.show()