import numpy as np

x = np.linspace(-6, 6, 2)
y = np.linspace(-6, 6, 2)

X, Y = np.meshgrid(x, y)
print(X)
print(Y)
