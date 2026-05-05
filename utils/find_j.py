# Authors: Konstantinos Papalamprou
# Date: 25/4/2026
#
# Introduction to ECE - team 09 - BJT
# Utility used ONCE to calculate jacobean matrix symbolically

# Libraries
import sympy as sp

# Symbols
V_1, V_2, V_3 = sp.symbols('V_1 V_2 V_3')
R_2 = sp.Symbol('R_2')

R_1, R_s, R_L = sp.symbols('R_1 R_s R_L')
I_s, V_t, V_a, V_cc = sp.symbols('I_s V_t V_a V_cc')
beta = sp.Symbol('beta')

# Symbolic functions
f1 = (V_1/R_2) + ((V_1 - V_cc)/R_1) + (I_s * (sp.exp((V_1 - V_2)/V_t) - 1))
f2 = (V_2/R_s) - (I_s * (sp.exp((V_1 - V_2)/V_t) - 1)) - (beta * I_s * (sp.exp((V_1 - V_2)/V_t) - 1) * (1 + ((V_3 - V_2) / V_a)))
f3 = ((V_3 - V_cc) / R_L) + (beta * I_s * (sp.exp((V_1 - V_2)/V_t) - 1) * (1 + ((V_3 - V_2) / V_a)))

F = sp.Matrix([f1, f2, f3])    # Matrix of functions
V = sp.Matrix([V_1, V_2, V_3]) # Matrix of Voltages

# Calculate the jcobian
J = F.jacobian(V)

print("J = ")
sp.pprint(J)

print("or in python syntax:")
for i in range(3):
    for j in range(3):
        print(f"df{j+1}_dV_{i+1} =", J[i,j])
    print("")