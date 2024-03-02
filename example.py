'''
This script provide an example to use PolyCtexpt.
'''


import numpy as np
import sympy as sp

from lib.PolyCtexpt import PolyCtexpt, poly_ctexpt_add, poly_ctexpt_matmul

Nt = 3
matrix_size = 2

# %% Define the expressions
rng = np.random.default_rng(0)
f = [PolyCtexpt(Nt+3) for kk in range(2)]
symb_f = [sp.zeros(matrix_size, matrix_size) for kk in range(2)]
symb_t = sp.symbols('t t1 t2', real=True)

# Generate the random coefficients
for nf in range(2):
    N_term = rng.integers(20, 100)  # How many terms in the polynomial?
    for n_term in range(N_term):
        # My implementation
        t_power = rng.integers(1, 5, Nt)
        t_expon = rng.uniform(-5, 5, Nt) + 1.j * rng.uniform(-10, 10, Nt)
        mat = rng.uniform(-100, 100, (matrix_size, matrix_size)) + 1.j * rng.uniform(-100, 100, (matrix_size, matrix_size))
        f[nf].insert_mat(mat, t_power, t_expon)

        # Sympy
        sp_power = 1
        sp_exp = 0
        for nt in range(Nt):
            sp_power *= symb_t[nt]**t_power[nt]
            sp_exp += sp.Float(np.real(t_expon[nt]))*symb_t[nt] + sp.I*sp.Float(np.imag(t_expon[nt]))*symb_t[nt] #t_expon[nt]*symb_t[nt]
        symb_f[nf] += sp.Matrix(mat) * sp_power * sp.exp(sp_exp)
        

# %% Test 1
t_subs = rng.uniform(-10, 10, Nt)
g_0 = poly_ctexpt_matmul(f[0].T(), f[1])
symb_g_0 = symb_f[0].transpose() @ symb_f[1]
g_0_result = g_0.evalf(t_subs)
symb_g_0_result = symb_g_0.evalf(subs={symb_t[kk]: t_subs[kk] for kk in range(Nt)})

# Relative error
print('Test 1')
print('Relative error of the matrix elements')
print('  ', np.abs(g_0_result-symb_g_0_result) / ((np.abs(g_0_result)**2 + np.abs(symb_g_0_result)**2)/2)**0.5)


# %% Test 2
t_subs = rng.uniform(-10, 10, Nt)
g_1 = poly_ctexpt_add(g_0, f[1])
symb_g_1 = symb_g_0 + symb_f[1]
g_1_result = g_1.evalf(t_subs)
symb_g_1_result = symb_g_1.evalf(subs={symb_t[kk]: t_subs[kk] for kk in range(Nt)})

# Relative error
print('Test 2')
print('Relative error of the matrix elements')
print('  ', np.abs(g_1_result-symb_g_1_result) / ((np.abs(g_1_result)**2 + np.abs(symb_g_1_result)**2)/2)**0.5)


# %% Test 3
t_subs = rng.uniform(-10, 10, Nt-1)
g_2 = g_1.substitute_t_1_2subt(from_idx=2, to_idx=(1,0))  # t_2 <- t_1 - t
g_2.simplify()
symb_g_2 = symb_g_1.subs(symb_t[2], symb_t[1]-symb_t[0])
g_2_result = g_2.evalf(t_subs)
symb_g_2_result = symb_g_2.evalf(subs={symb_t[kk]: t_subs[kk] for kk in range(Nt-1)})

# Relative error
print('Test 3')
print('Relative error of the matrix elements')
print('  ', np.abs(g_2_result-symb_g_2_result) / ((np.abs(g_2_result)**2 + np.abs(symb_g_2_result)**2)/2)**0.5)
