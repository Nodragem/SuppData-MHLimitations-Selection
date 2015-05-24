import matplotlib.pyplot as plt
import numpy as np


def MH(A, mu, sigma, K, B, x):
    g1 = (1+B)*np.exp(-((x-mu)**2)/(2*(sigma)**2)) #+ A*exp((-(x-mu-100)**2)/(sigma**2))
    g2 = B*np.exp(-((x-mu)**2)/(2*(K*sigma)**2))
    return A*(g1-g2)

x = np.arange(0, 30, 0.1)
SA1 = MH(200, 0, 5.0, 1.2, 6.0,x)
SA2 = MH(200, 0, 5.0, 2.0, 1.42868,x)
SA3 = MH(200, 0, 5.0, 1.2, 8.0,x)

plt.plot(x, SA1, "-", x, SA2, "r--", x, SA3, "g-.")
plt.show()