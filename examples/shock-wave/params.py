Mach = 2
g = 5. / 3.
Ru = 8.3144598 # Universal gas constant
Mol = 40e-3
Rg = Ru / Mol 

n_in=1.6095e+21
T_in=300.0
u_in = Mach * ((g * Rg * T_in) ** 0.5)

n_out = (g + 1.) * Mach * Mach / ((g - 1.) * Mach * Mach + 2.) * n_in
u_out = ((g - 1.) * Mach * Mach + 2.) / ((g + 1.) * Mach * Mach) * u_in
T_out = (2. * g * Mach * Mach - (g - 1.)) * ((g - 1.) * Mach * Mach + 2.) / ((g + 1) ** 2 * Mach * Mach) * T_in

print('u_in =', u_in)

print('n_out =', n_out)
print('u_out =', u_out)
print('T_out =', T_out)
