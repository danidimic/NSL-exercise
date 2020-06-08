import math
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

a = 1.5
sigma = 0.20 #standard deviation of noise
Nepochs = 15 #numero di epoche
Ntrain = 1000 #numero di valori di test
Nvalid = 100  #numero di valori di validazione

name = ''

def target(X):
	x = X[0]
	y = X[1]
	return np.sin(x**2 + y**2)

#generate training inputs
np.random.seed(0)

# generate synthetic training dataset
x_train = np.zeros((Ntrain, 2))
y_train = np.zeros(Ntrain)
for i in range(Ntrain):
	x_train[i,0] = np.random.uniform(-a, a)
	x_train[i,1] = np.random.uniform(-a, a)
	y_train[i] = target(x_train[i])


# generate synthetic validation data
x_valid = np.zeros((Nvalid, 2))
y_valid = np.zeros(Nvalid)
for i in range(Nvalid):
	x_valid[i,0] = np.random.uniform(-a, a)
	x_valid[i,1] = np.random.uniform(-a, a)
	y_valid[i] = target(x_valid[i])


#compose the NN model
model = Sequential()

model.add(Dense(units=27, input_dim=2, activation='relu'))
model.add(Dense(units=18, activation='relu'))
model.add(Dense(10, activation='softmax'))

#compile the model choosing optimizer, loss and metrics objects
model.compile(optimizer='sgd', loss='mse', metrics=['mse'])
#get a summary of our composed model
model.summary()

# fit the model using training dataset
history = model.fit(x=x_train, y=y_train, 
          batch_size=32, epochs=Nepochs,
          shuffle=True, # a good idea is to shuffle input before at each epoch
          validation_data=(x_valid, y_valid))

'''
x_predicted = np.linspace(-1.5, 1.5, 100)
y_predicted = model.predict(x_predicted).reshape(-1)

plt.scatter(x_predicted, y_predicted, color='tab:orange', label='prediction') #predizione NN
plt.scatter(x_valid, y_valid, color='r', label='validation data') #dati per validazione
plt.plot(xtarg, target(xtarg), label='target', lw=2.5) #funzione target
plt.grid(True)
plt.legend()
plt.show()
'''
n = 10
x_predicted = np.zeros((n, 2))
for i in range(n):
	x_predicted[i,0] = np.random.uniform(-a, a)
	x_predicted[i,1] = np.random.uniform(-a, a)
y_predicted = model.predict(x_predicted).reshape(-1)

print(x_predicted)
print()
print()
print(y_predicted)
print(y_predicted.size())

x = x_predicted[:,0]
y = x_predicted[:,1]
z = y_predicted[:]

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, c='tab:blue')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


#save_model_path='exercise3/model' + name
#model.save(filepath=save_model_path, include_optimizer=True)


