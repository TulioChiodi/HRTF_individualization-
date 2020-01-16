# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:20:34 2020

@author: rdavi
"""

from joblib import Parallel, delayed
import numpy as np
import scipy.io as sio

#-- KERAS
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras import optimizers
from keras.callbacks import EarlyStopping

#%% DATALOADER
# Data from matlab script DeepNet.m
path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas\Python\data_for_train.mat'
mat_contents = sio.loadmat(path)

inpt = mat_contents['input']
target = mat_contents['target']
sig = mat_contents['sig']
mu = mat_contents['mu']

(no_azi, no_ele, no_samles, no_subjects, no_channels) = np.shape(target)

#%% BUILD NETWORK
def network(x, y):
    ### LAYERS ###
    nodes = 32
    model = Sequential()
    
    model.add(Dense(nodes, activation='linear', input_dim = 37,
                     kernel_initializer='glorot_uniform',
                     bias_initializer='zeros'))
    model.add(Dropout(0.5))
    
    model.add(Dense(nodes, activation='relu',
                     kernel_initializer='glorot_uniform',
                     bias_initializer='zeros'))
    model.add(Dropout(0.5))
    
    model.add(Dense(nodes, activation='relu',
                     kernel_initializer='glorot_uniform',
                     bias_initializer='zeros'))
    model.add(Dropout(0.5))
    
    model.add(Dense(nodes, activation='relu',
                     kernel_initializer='glorot_uniform',
                     bias_initializer='zeros'))
    model.add(Dropout(0.5))
    
    model.add(Dense(nodes, activation='relu',
                     kernel_initializer='glorot_uniform',
                     bias_initializer='zeros'))
    model.add(Dropout(0.5))
        
    model.add(Dense(200, activation='relu'))
       
    ### Compile (OPTIONS) ###
    adam = optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, amsgrad=False)
    model.compile(optimizer = adam,
                  loss = 'mse')
    es = EarlyStopping(monitor='val_loss', mode='min', verbose=2, patience=500)
    ### Training ###
    model.fit(x, y, 
              epochs=20000, 
              verbose=0,            
              validation_split=0.1,  
              shuffle=True, 
              initial_epoch = 0,
              workers=1,
              callbacks=[es])
    return model

#%% Training 
# MULTIPROCESSING thread
def elev(l,k):
    for m in range(0, no_ele):
        print("Azimuth:"+str(l)+" Elevação:"+str(m)+" Canal:"+str(k))
        xi = np.transpose(inpt[:,:,k])
        yi = np.transpose(target[l,k,:,:,k])
        net = network(xi, yi)
        #Save network
        net.save("Deep_Keras_C_64/net_"+str(l+1)+"_"+str(m+1)+"_"+str(k+1)+".h5")

# Parallel azyymuth
for k in range(0, no_channels):
    Parallel(n_jobs=2, prefer="threads")(
            delayed(elev(l, k)) for l in range(0, no_azi))


    
















