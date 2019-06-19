import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import pandas as pd
from scipy import interpolate

x = [4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
y = [2.63, 2.14, 1.77, 1.44, 1.23, 1.06, 0.941, 0.38, 0.21, 0.13, 0.087, 0.063, 0.048, 0.038, 0.03, 0.025]

f = interpolate.interp1d(x,y,kind='cubic')

xnew = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,45,55,65,75,85,95,100]
ynew = f(xnew)
for i in range(12):
	print(xnew[i], ynew[i])

plt.plot(x,y, 'o', xnew, ynew, '-')
plt.show()
