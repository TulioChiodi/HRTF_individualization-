# HRTF_individualization

Preprocess_CIPIC_ARI_ITA.m - The main preprocessing code, it makes the assembly of CIPIC, ARI and ITA HRIR datasets.
                             As these dataasets where measured using different setups, here we put them in the 
                             same coordinate system, sample rate, IR length, find the subjects with complete anthropometry,
                             transform from HRIR -> HRTF -> DTF, filtering the DTF between 20 - 18kHz. 
                             The product of this code is a 5-D matrix containing DTFs for 118 subjects, for both ears, 25 azymuth and 50                                elevations (mainly the CIPIC data structre).
                             
                               
                           
