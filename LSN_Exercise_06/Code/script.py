import os

nvalues = 100
add = 1.5/nvalues

for i in range(nvalues):
	print("blocco: ", i+1)
	#esegui programma c++
	cmd = './Monte_Carlo_ISING_1D.exe'
	os.system(cmd)
	#leggi temperatura in input.dat
	lines = open("input.dat", 'r').readlines()
	t = float(lines[0])
	t += add
	lines[0] = str(t)
	lines[0] += "\n"
	#modifica temperatura in input.dat
	out = open("input.dat", 'w')
	out.writelines(lines)
	out.close()
