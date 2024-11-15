import time
import numpy as np


beta0 = np.deg2rad(90.0)       # beta0  = 90 deg
l1 = 10
start = time.time()
for i in range(10000000):
    a = l1 * np.exp(1j * beta0)
end = time.time()
print(end-start)