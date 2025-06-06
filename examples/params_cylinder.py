from math import pi as PI

Mach = 10.0
g = 5. / 3.
Ru = 8.3144598 # Universal gas constant
Mol = 40e-3
Rg = Ru / Mol 

Na = 6.02214129e+23
m = Mol / Na

u_in = 2624.0
T_in = (u_in ** 2) / (Mach**2 * g * Rg)
rho_in = 1.127e-6
n_in = rho_in / m
delta = 1.6
Kn = 8. / (5 * PI**.5 * delta)

print('n =', n_in)
print('u =', u_in)
print('T =', T_in)
print('Kn=', Kn)
