# HRTF_individualization

Preprocess_CIPIC_ARI_ITA.m - The main preprocessing code, it makes the assembly of CIPIC, ARI and ITA HRIR datasets.
                             As these dataasets where measured using different setups, here we put them in the 
                             same coordinate system, sample rate, IR length, find the subjects with complete anthropometry,
                             transform from HRIR -> HRTF -> DTF, filtering the DTF between 20 - 18kHz. 
                             The product of this code is a 4-D matrix containing DTFs for 118 subjects, for both ears, 25     
                             azymuth and 50 elevations (mainly the CIPIC data structre).
                             
                               
PCA_CIPIC_ARI_ITA.m        - Performs principal component analysis over the structured data from 'Preprocess_CIPIC_ARI_ITA.m'



NeuralNet_CIPIC.m          - Train shallow neural network with anthropometric data as input features and principal components as targets.                              (for the CIPIC dataset only).


NeuralNet_CIPIC_ARI_ITA.m  - Train shallow neural network with anthropometric data as input features and principal components as targets.                              (for the CIPIC, ARI and ITA dataset together).

