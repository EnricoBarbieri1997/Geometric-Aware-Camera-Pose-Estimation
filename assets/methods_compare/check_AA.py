import sympy as sp

# Define symbols
fx, fy, cx, cy = sp.symbols('f_x f_y c_x c_y', real=True)

# Intrinsic matrix with zero skew
A = sp.Matrix([
    [fx,  0,  cx],
    [0,  fy,  cy],
    [0,   0,   1]
])

# Compute omega = A^{-T} * A^{-1}
A_inv = A.inv()
omega = A_inv.T * A_inv

# Simplify entries
omega = sp.simplify(omega)

# Display result
print("ω = A^{-T} * A^{-1}:")
sp.pprint(omega)
