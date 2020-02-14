# Synthesis of individualized HRTF with neural networks

The main idea of the project is to synthesise individualized HRTFs from anthropometric measurements, with the use initially of shallow neural networks, and later experiments with deeper networks will be performed.
Here we use three main databases to expand training data, the CIPIC, ARI and ITA databases (http://sofacoustics.org/data/database/).


  ## Main codes 
### Preprocess_CIPIC_ARI_ITA.m   

The main preprocessing code, it makes the assembly of CIPIC, ARI and ITA HRIR datasets.
As these dataasets where measured using different setups, here we put them in the 
same coordinate system, sample rate, IR length, find the subjects with complete anthropometry,
transform from HRIR -> HRTF -> DTF, filtering the DTF between 20 - 18kHz. 
The product of this code is a 4-D matrix containing DTFs for 118 subjects, for both ears, 25     
azymuth and 50 elevations (mainly the CIPIC data structre).
                     
                               
### PCA_CIPIC_ARI_ITA.m   
 
Performs principal component analysis over the structured data from 'Preprocess_CIPIC_ARI_ITA.m'
  


### NeuralNet_CIPIC.m    

Train shallow neural network with anthropometric data as input features and principal components as targets.(for the CIPIC dataset only).
 

### NeuralNet_CIPIC_ARI_ITA.m     

Train shallow neural network with anthropometric data as input features and principal components as targets.(for the CIPIC, ARI and ITA dataset together).


### Anthropometry_CIPIC_ARI_ITA.m  

Assemble the anthropometry features to train CIPIC_ARI_ITA networks.


## Built With
* [Matlab](https://www.mathworks.com/products/matlab.html)
