import os
import numpy as np

nvalues = 100
temp = np.linspace(0.5, 2.0, nvalues)

for i in range(nvalues):

	print("blocco: ", i+1)
	#leggi temperatura in input.dat
	lines = open("input.dat", 'r').readlines()
	lines[0] = str(temp[i])
	lines[0] += "\n"
	#modifica temperatura in input.dat
	out = open("input.dat", 'w')
	out.writelines(lines)
	out.close()
	#esegui programma c++
	cmd = './Monte_Carlo_ISING_1D.exe'
	os.system(cmd)
