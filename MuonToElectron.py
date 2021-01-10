
import numpy as np

x1 = 1
third = 1 / 3
f = 4 / np.pi # Ettore Majorana's ring of a disk potential factor.
f2 = f * f

# Iterate to most stable root.
for i in range(2000):
    x1 = np.power((192 * x1 * x1 + 2 * x1 * f - f2) / 192, third)

a = 1/(x1 - 1) # Negative charge.

print('Xi = 4/Pi, a = %.48f' % a)

x3 = 1
x4 = 1
f = 95 / 96
f2 = f * f

# Iterate to most stable roots.
for i in range(2000):
    x3 = np.power((192 * x3 * x3 + 2 * x3 * f - f2) / 192, third)
    x4 = np.power((192 * x4 * x4 - 2 * x4 * f - f2) / 192, third)

c = 1/(x3 - 1) # Negative charge.
d = 1/(1 - x4) # Positive charge.

print('Xi = 95/96, c = %.48f, d = %.48f' % (c, d))
print('Xi = 95/96, c * d = %.48f' % (c * d))

print('Approximated mass ratio between the Muon and the electron %.48f'
      % (a * (1 + (x3-1)*(1-x4))))
'''
Output:
Xi = 4/Pi, a = 206.751339885031512721980107016861438751220703125000
Xi = 95/96, c = 192.046394360137696821766439825296401977539062500000, d = 63.541358229208128705067792907357215881347656250000
Xi = 95/96, c * d = 12202.888740665284785791300237178802490234375000000000
Approximated mass ratio between the Muon and the electron 206.768282704414616546273464336991310119628906250000
'''
