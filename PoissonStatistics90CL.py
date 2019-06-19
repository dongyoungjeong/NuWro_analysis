import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy import stats
from math import factorial

N_b_el = [[1.0076727, 0.8182366, 0.5985544, 0.5619407, 0.4823457, 0.382056, 0.334299, 0.3008691, 0.2833582, 0.2817663, 0.2133146, 0.2101308], # sorted by energy
                [0.9981213, 0.8293799, 0.6351681, 0.6096977, 0.5269189, 0.4345887, 0.3947912, 0.3374828, 0.31838, 0.3120124, 0.2881339, 0.2610716],
                [1.0952272, 0.8930559, 0.6590466, 0.668598, 0.5714921, 0.5046323, 0.4648348, 0.3947912, 0.3740965, 0.3215638, 0.3295233, 0.3120124]]
N_b_mu = [[3.7428679, 3.5949934, 3.2795278, 2.9607761, 2.7866128, 2.5204387, 2.2641229, 2.1096762, 1.9782322, 1.8204994, 1.6989137, 1.6069029], # sorted by energy
            [3.7724428, 3.516127, 3.3321054, 3.0199259, 2.7997572, 2.5730163, 2.4218557, 2.2378341, 2.0603847, 1.8927936, 1.7974967, 1.6956276],
            [3.7231513, 3.5917073, 3.4569772, 3.0725035, 2.9114846, 2.6978881, 2.4744333, 2.3298449, 2.1261067, 1.9618017, 1.8435021, 1.7482052]]
N_b_combined = [[4.321908, 3.921912, 3.351186, 3.058506, 2.807289, 2.456073, 2.192661, 2.026809, 1.90242, 1.782909, 1.587789, 1.514619], # sorted by energy
                [4.329225, 3.880449, 3.446307, 3.175578, 2.885337, 2.575584, 2.402415, 2.178027, 2.017053, 1.882908, 1.775592, 1.65852],
                [4.441419, 4.034106, 3.575574, 3.304845, 3.036555, 2.775582, 2.548755, 2.334123, 2.151198, 1.948761, 1.873152, 1.775592]]
eff_el = [[0.063336, 0.095787, 0.128712, 0.15978, 0.18856, 0.212684, 0.234872, 0.252656, 0.269045, 0.281724, 0.29235, 0.302056],
                 [0.06191, 0.093787, 0.124296, 0.154956, 0.182336, 0.20628, 0.227471, 0.24586, 0.262045, 0.275724, 0.28674, 0.297056],
                 [0.059058, 0.088524, 0.115632, 0.142308, 0.166112, 0.18726, 0.205867, 0.222268, 0.23607, 0.247152, 0.25757, 0.26628]]
eff_mu = [[0.058206, 0.087682, 0.117216, 0.144896, 0.170704, 0.193856, 0.213664, 0.231472, 0.24528, 0.258366, 0.2684, 0.278116],
                 [0.05678, 0.085261, 0.1128, 0.13966, 0.163888, 0.186452, 0.205862, 0.223278, 0.237095, 0.24958, 0.26023, 0.27034],
                 [0.053354, 0.079419, 0.104552, 0.128424, 0.150848, 0.170028, 0.187258, 0.202084, 0.214515, 0.225008, 0.23467, 0.243176]]

def poisson(k,lam):
    return lam**k*np.exp(-lam)/factorial(k)
def N_s_cl90(N_b):
    x = np.linspace(0, 100, 100000)
    y = np.array([poisson(round(N_b), xx) for xx in x])
    diff_min = 1.0
    cl90_index = 0
    for i in range(0, len(x)):
        diff = abs(y.cumsum()[i]/y.cumsum().max() - 0.9)
        if diff < diff_min:
            diff_min = diff
            cl90_index = i
    return x[cl90_index]-N_b
def printN_s(N_b_list):
    N_s = np.array([lim_90(xx) for xx in N_b_list])
    print(N_s)

def lim_90(k):
    return scipy.optimize.brentq(lambda x: scipy.stats.poisson.cdf(round(k), x) - .1, 0, 2*k+3) - k
def combined_90(r, c):
    return scipy.optimize.brentq(lambda x: scipy.stats.poisson.cdf(round(N_b_el[r][c]), N_b_el[r][c] + eff_el[r][c] * x) * scipy.stats.poisson.cdf(round(N_b_mu[r][c]), N_b_mu[r][c] + eff_mu[r][c] * x) - .1, 0, 100)

#print("Electron-neutrino")
#for i in range(0,12):
#    printN_s(N_b_electron[i])
print("Muon-neutrino")
for i in range(len(N_b_mu)):
    printN_s(N_b_mu[i])
print("Electron & Muon combined")
for i in range(3):
    a = np.array([combined_90(i, xx) for xx in range(12)])
    print(list(a))