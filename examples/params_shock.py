from math import pi as PI

Mach = 8
g = 5. / 3.
Ru = 8.3144598 # Universal gas constant
Mol = 40e-3
Rg = Ru / Mol 

Na = 6.02214129e+23
m = Mol / Na

n_in=2e+23
T_in=200.0
u_in = Mach * ((g * Rg * T_in) ** 0.5) # was this
# u_in = Mach * ((g * T_in) ** 0.5)

n_out = (g + 1.) * Mach * Mach / ((g - 1.) * Mach * Mach + 2.) * n_in
u_out = ((g - 1.) * Mach * Mach + 2.) / ((g + 1.) * Mach * Mach) * u_in
T_out = (2. * g * Mach * Mach - (g - 1.)) * ((g - 1.) * Mach * Mach + 2.) / ((g + 1) ** 2 * Mach * Mach) * T_in

Kn = 0.564
delta = 8. / (5 * PI**.5 * Kn)

mu_0 = 2.125e-5
T_0  = 273.11
C    = 144.4

mu_suth = lambda T: mu_0 * ((T_0 + C) / (T + C)) * ((T / T_0) ** (3. / 2.))

l_s = (delta * mu_suth(T_in) * (2. * Rg * T_in)**0.5) / (m * n_in * Rg * T_in)

print('l_s =', l_s)

print('n_in  =', n_in)
print('u_in  =', u_in)
print('T_in  =', T_in)

print('n_out =', n_out)
print('u_out =', u_out)
print('T_out =', T_out)

print('Mach_paper =', u_in / (g * T_in) ** .5)
