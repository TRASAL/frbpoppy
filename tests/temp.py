import math

f = 1.4e9
bw = 0.1e9
f_1 = f - 0.5*bw
f_2 = f + 0.5*bw
si = 0
sp = si + 1
sm = si - 1
z = 0.75
lum = 8e44 * 1e-7
pc = 3.086e16
dist = 2.68214560e9*pc
f = 1.4 * 1e6
f_low = 10e6
f_high = 10e9

freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)
nom = lum * (1+z)**sm * freq_frac
den = 4*math.pi*dist**2 * (f_high**sp - f_low**sp)
s_peak = nom/den
s_peak *= 1e26
print(s_peak)
