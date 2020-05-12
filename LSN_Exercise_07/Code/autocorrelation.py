import numpy as np
import statsmodels.tsa.stattools as stat

M = 5e5
n = 100

phase = 'gas'
ene = np.loadtxt(phase + '_ene.out')
pres = np.loadtxt(phase + '_pres.out')

Aene = stat.acf(ene, fft=True, nlags=n-1)
Apres = stat.acf(pres, fft=True, nlags=n-1)

np.savetxt('Autocorr_ene.out', Aene)
np.savetxt('Autocorr_pres.out', Apres)
