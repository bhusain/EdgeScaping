import numpy as np
import matplotlib.pyplot as plt

from keras.models import Sequential
from keras.layers.core import Dense
from keras.optimizers import SGD

from sklearn.model_selection import train_test_split

from time import time
import numpy as np
import keras.backend as K
from keras.engine.topology import Layer, InputSpec
from keras.layers import Dense, Input
from keras.models import Model
from keras.optimizers import SGD
from keras.utils.vis_utils import plot_model

nb_classes = 10


text_file = open("output_new.txt", "r")
lines = text_file.readlines()
#print lines
print(len(lines))
text_file.close()

data = []
for i in lines:
    i = i.split('[', 1)[-1]
    i = i[:-2]
    line = np.fromstring(i, dtype=int, sep=' ')
    data.append(line)

print(len(data))

data_train = data

#data_train = data[0:3500]
#data_test = data[3501:]


model = Sequential()
model.add(Dense(512, activation='relu', input_shape=(784,)))
model.add(Dense(512, activation='relu'))
model.add(Dense(10, activation='softmax'))

model.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.001),
              metrics=['accuracy'])

data_train, data_val= train_test_split(data_train)
dims=[361, 500, 500, 2000, 10]

#print(data_train)

print(np.array(data_train).shape)

def autoencoder(dims, act='relu'):
    """
        Fully connected auto-encoder model, symmetric.
        Arguments:
        dims: list of number of units in each layer of encoder. dims[0] is input dim, dims[-1] is units in hidden layer.
        The decoder is symmetric with encoder. So number of layers of the auto-encoder is 2*len(dims)-1
        act: activation, not applied to Input, Hidden and Output layers
        return:
        Model of autoencoder
        """
    n_stacks = len(dims) - 1
    # input
    x = Input(shape=(dims[0],), name='input')
    h = x
    
    # internal layers in encoder
    for i in range(n_stacks-1):
        h = Dense(dims[i + 1], activation=act, name='encoder_%d' % i)(h)

    # hidden layer
    h = Dense(dims[-1], name='encoder_%d' % (n_stacks - 1))(h)  # hidden layer, features are extracted from here

    # internal layers in decoder
    for i in range(n_stacks-1, 0, -1):
        h = Dense(dims[i], activation=act, name='decoder_%d' % i)(h)
    
    # output
    h = Dense(dims[0], name='decoder_0')(h)
    
    return Model(inputs=x, outputs=h)

model = autoencoder(dims)
model.compile(loss='kld', optimizer=SGD(lr=0.01, momentum=0.9))

print(model.summary())
model.fit(np.array(data_train), batch_size=20, epochs=2, verbose=1)






