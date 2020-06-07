import math
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

m = 2  #slope
b = 1  #intersect

sigma = 0.25 #standard deviation of noise
Nepochs = 100 #numero di epoche
Ntrain = 1000 #numero di valori di test

def target(x):
	return m*x + b

#generate training inputs
np.random.seed(0)
x_valid = np.random.uniform(-1, 1, 10)
x_valid.sort()
y_target = target(x_valid)
x_train = np.random.uniform(-1, 1, Ntrain)

#random noise
y_train = np.random.normal( target(x_train), sigma ) # actual measures from which we want to guess regression parameters
y_valid = np.random.normal( target(x_valid), sigma )

#compose the NN model
model = tf.keras.Sequential()
model.add(Dense(1, input_shape=(1,)))

#compile the model choosing optimizer, loss and metrics objects
model.compile(optimizer='sgd', loss='mse', metrics=['mse'])
#get a summary of our composed model
model.summary()

# fit the model using training dataset
history = model.fit(x=x_train, y=y_train, 
          batch_size=32, epochs=Nepochs,
          shuffle=True, # a good idea is to shuffle input before at each epoch
          validation_data=(x_valid, y_valid))


# return weights
model.get_weights()

# evaluate model
score = model.evaluate(x_valid, y_valid, batch_size=32, verbose=1)
# print performance
print()
print('Test loss:', score[0])
print('Test accuracy:', score[1])

# evaluate model with the exact curve
score = model.evaluate(x_valid, y_target, batch_size=32, verbose=1)
# print performance
print()
print('Test loss:', score[0])
print('Test accuracy:', score[1])

x_predicted = np.random.uniform(-1, 1, 100)
y_predicted = model.predict(x_predicted)

'''
plt.title("Dati di validazione")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x_valid, y_target)
plt.scatter(x_valid, y_valid, color='r')
plt.grid(True)
plt.show()

plt.title("Regressione rete neurale")
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(x_predicted, y_predicted, color='r', label='predizione')
plt.plot(x_valid, y_target, label='target')
plt.legend()
plt.grid(True)
plt.show()
'''

y_pred = []
for i in range (y_predicted.size):
	y_pred.append(y_predicted[i][0])

#out = 'exercise1/Nepochs' + str(Nepochs) + '.out'
#out = 'exercise1/Ntrain' + str(Ntrain) + '.out'
out = 'exercise1/sigma' + str(sigma) + '.out'
outval = 'exercise1/valid' + str(sigma) + '.out'
np.savetxt(out, np.column_stack([x_predicted, y_pred]) )
np.savetxt(outval, np.column_stack([x_valid, y_valid]) )

