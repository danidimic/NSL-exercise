import os
import math
import numpy as np

#cancella risultati precedenti
cmd = './clean.sh'
os.system(cmd)

nvalues = 200
mu = np.linspace(0.8, 0.9, nvalues)
sigma = np.linspace(0.6, 0.7, nvalues)

for i in range(nvalues):

	print("ciclo: ", i+1)
	print()
	#leggi valori in input.dat
	lines = open("input.dat", 'r').readlines()
	lines[1] = str(mu[i])
	lines[1] += "\n"

	for j in range(nvalues):
		lines[2] = str(sigma[j])
		lines[2] += "\n"
		#modifica sigma, mu in input.dat
		out = open("input.dat", 'w')
		out.writelines(lines)
		out.close()
		#esegui programma c++
		cmd = './exercise_1.exe'
		os.system(cmd)


mu, sigma, H = np.loadtxt("results/minimize.out", unpack=True)

index = np.argmin(H)

print("Valori mu, sigma che minimizzano H")
print("mu = ", mu[index])
print("sigma = ", sigma[index])

print()
print("Energia del ground state", H[index])
