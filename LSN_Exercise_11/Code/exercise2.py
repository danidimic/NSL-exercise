import math
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

sigma = 0.20 #standard deviation of noise
Nepochs = 150 #numero di epoche
Ntrain = 1000 #numero di valori di test

name = ''

def target(x):
	return 4 - 3*x - 2*x**2 + 3*x**3

#generate training inputs
np.random.seed(0)

#validation data
x_valid = np.linspace(-1, 1, 15)
x_valid.sort()
y_valid = np.random.normal( target(x_valid), sigma )
#train data
x_train = np.random.uniform(-1, 1, Ntrain)
y_train = np.random.normal( target(x_train), sigma ) # actual measures from which we want to guess regression parameters
#funzione target
xtarg = np.linspace(-1,1, 25)
y_target = target(xtarg)

#compose the NN model
model = Sequential()

model.add(Dense(units=200, input_dim=1, activation='relu'))
model.add(Dense(units=45, activation='relu'))
model.add(Dense(units=1))

#compile the model choosing optimizer, loss and metrics objects
model.compile(optimizer='sgd', loss='mse', metrics=['mse'])
#get a summary of our composed model
model.summary()

# fit the model using training dataset
history = model.fit(x=x_train, y=y_train, 
          batch_size=32, epochs=Nepochs,
          shuffle=True, # a good idea is to shuffle input before at each epoch
          validation_data=(x_valid, y_valid))


x_predicted = np.linspace(-1.5, 1.5, 100)
y_predicted = model.predict(x_predicted).reshape(-1)

plt.scatter(x_predicted, y_predicted, color='tab:orange', label='prediction') #predizione NN
plt.scatter(x_valid, y_valid, color='r', label='validation data') #dati per validazione
plt.plot(xtarg, target(xtarg), label='target', lw=2.5) #funzione target
plt.grid(True)
plt.legend()
plt.show()

save_model_path='exercise2/model' + name
model.save(filepath=save_model_path, include_optimizer=True)





