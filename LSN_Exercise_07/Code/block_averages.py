import math
import numpy as np
import matplotlib.pyplot as plt


def data_blocking(arr, N, L):
	ave = []
	err = []
	#Determino i valori medi su N blocchi di lunghezza L
	for i in range(N):
		s = 0
		for j in range(L):
			s  += arr[i*L+j]
		s  /= L
		ave.append(s)

	#Determino la media e l'incertezza su tutti i blocchi
	s = 0
	s2 = 0
	for i in range(N):
		s += ave[i]
		s2 += ave[i]**2		

	s  /= N
	s2 /= N
	
	average  = s
	error = math.sqrt( (s2-s**2)/(N-1) )
	return [average, error]

phase = 'gas'
ene  = np.loadtxt(phase + "_ene.out")
pres = np.loadtxt(phase + "_pres.out")

aveE = []
errE = []
aveP = []
errP = []

M = 5e5 #valori totali
L = np.arange(start=10, stop=5000, step=25)

#L1 = np.arange(start=10, stop=100, step=10)
#L2 = np.arange(start=100, stop=1000, step=20)
#L3 = np.arange(start=1000, stop=5000, step=200)

#L12 = np.concatenate( (L1,L2) )
#L = np.concatenate( (L12,L3) )
#L = np.append( L, 5000 )

for i in range(L.size):
	print("ciclo :", i)
	N = int(M/L[i]) #numero di blocchi

	if(N*L[i]>M):
		print(i)

	db = data_blocking(ene, N, L[i])
	aveE.append(db[0])
	errE.append(db[1])

	db = data_blocking(pres, N, L[i])
	aveP.append(db[0])
	errP.append(db[1])	

d = np.column_stack((L, aveE))
E =  np.column_stack((d, errE))
d = np.column_stack((L, aveP))
P =  np.column_stack((d, errP))

np.savetxt('DB_ene.out', E, delimiter=" " )
np.savetxt('DB_pres.out', P, delimiter=" " )
