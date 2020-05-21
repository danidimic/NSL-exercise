import os
import math
import numpy as np


enemin = 10
muopt = 10
sigmaopt = 10

ntot = 0.
naccept = 0.
nstep = 25
delta = 0.015

#cancella risultati precedenti
cmd = './clean.sh'
os.system(cmd)

#temperatura
t1 = np.arange(1, 15.5, step=0.5)
t2 = np.arange(0.1, 1, step=0.1)
temp = np.concatenate((t2,t1))
temp = np.flip(temp)

beta = 1./temp
current_pos = np.array([0.8, 0.6])
proposed_pos = np.array([])


#Calcolo variazionale tramite codice c++
def variationalMC(pos):
	lines = []
	#valori mu,sigma per calcolo
	lines.append(str(2.55) + "\n")
	lines.append(str(pos[0]) + "\n")
	lines.append(str(pos[1]) + "\n")
	#modifica sigma, mu in input.dat
	out = open("input.dat", 'w')
	out.writelines(lines)
	out.close()
	#esegui programma c++
	cmd = './variationalMC.exe'
	os.system(cmd)
	#leggi energia corrispondente a mu,sigma
	ene = np.loadtxt("results/minimize.out")
	return ene


#Ciclo per minimizzare mu, sigma
for i in range( len(beta) ):
	print("ciclo: ", i+1, "/", len(beta))
	print("temperatura fittizia", 1./beta[i])
	print()

	for j in range(nstep):

		#posizione proposta del Metropolis
		transition = np.random.uniform(low=delta, high=delta, size=2)
		proposed_pos = current_pos + transition;

		#energia della posizione corrente
		current_ene = variationalMC(current_pos)
		#energia della posizione proposta
		proposed_ene = variationalMC(proposed_pos)

		deltaE = proposed_ene - current_ene
		alfa = min( 1, np.exp(-beta[i]*deltaE) )
		r = np.random.uniform()

		if r < alfa:
			current_pos = proposed_pos	
			naccept += 1
		ntot += 1

		if current_ene < enemin:
			enemin = current_ene
			muopt = current_pos[0]
			sigmaopt = current_pos[1]

print()
print("Accettazione Metropolis: ", naccept/ntot )

print()
print("Energia del ground state", enemin)

print()
print("parametri variazionali ottimali:")
print("mu = ", muopt)
print("sigma = ", sigmaopt)
