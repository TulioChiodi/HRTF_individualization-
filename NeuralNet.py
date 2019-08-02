# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 04:07:30 2019

@author: rdavi
"""

# Convert matlab input nd target matrix to python
import scipy.io as sio
import numpy 

path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas\DADOS_TREINAMENTO\python_convert.mat'
mat_contents = sio.loadmat(path)
weights = mat_contents['weights']
InptMtx = mat_contents['InputMatrix']

(no_PC, no_subjects, no_directions, no_channels) = numpy.shape(weights)