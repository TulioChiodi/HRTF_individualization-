# DOCUMENTAÇÃO PT-BR
# Funções auxiliares à sinteze de HRTF individualizadas
Maior parte das funções são de autoria própria, mas algumas delas são de diversos autores, 
os devidos créditos são dados em cada algoritmo. 

Algumas da rotinas não estão sendo mais utilizadas nos códigos principais por terem sido encontrados melhores abordagens para o problema, mas foram de fundamental importancia para o desenvolvimento como um todo.


 ## Funções
### AA_ARI2CIPIC.m 
Realiza a transformação de hrtfs do banco de dados ARI em seu formato de padronização original para o formato de padronização adotado no banco CIPIC (ambos atualmente obsoletos pela adoção do formato SOFA). 


### hor2geo.m
Conversão de coordenadas horizontais polares para geodésicas.


### hor2sph.m
Conversão de coordenadas horizontais polares para esférica.


### hrtf2DTF.m
Transformação de hrtf em dtf pelo método de Middlebrooks (média logaritmica), Esta função é semelhante à SOFAhrtf2dtf() encontrada no SOFA-Toolbox, recomenda-se a utilização da função do Toolbox por possuir menos bugs. 


### itd_metodo3.m
Sintese itd a partir da largura e profundidade da cabeça (modelo Woodsworth otimizado por Algazi), aproximação esférica. 


### itd_metodo4.m
Calculo do itd a partir das relações de fase entre as HRTFs de ambas as orelhas (útil para HRTFs medidas, pois possuem naturalmente informação de fase).


### LSD.m 
Calculo da distorção espectral logaritmica entre duas hrtfs (recomenda-se o uso da função spec_dist.m, por ser um pouco mais robusta)


### minphase.m 
Calcula a fase mínima a partir da magnitude apenas e aplica a fase de excesso (ITD) sobre a RI [ITA-Toolbox] (o uso do Toolbox torna o processo bem mais lento, recomenda-se o uso da função phase_job.m para realizar a mesma operação).	


### nav2sph.m
Conversão de coordenadas navegacionais para coordenadas esfericas.


### nn_input_preprocess.m 
Normalização das linhas de uma matriz para terem média zero e variancia um (alternativa aao uso das funções built-in do matlab mapstd() ou zscore()).
Esse processo é útil para que as features estejam numa mesma escala e assim evitar a saturação nas funções de ativação da rede e para que os gradientes fluam no processo de backprop.


### PCA_fun.m
Realiza a analise de componentes principais pelo metodo de decomposição 
dos autovalores a partir da matriz de covariancia. Corresponde ao uso da função built-in do matlab pca() no modo 'eig'. 


### phase_job.m 
Calculo da fase mínima a partir de magnitude apenas e implemtação de atraso temporal correspondente ao ITD.


### pT_nOctaveBands.m 
Subdivisão de espectro em n-bandas de oitava. 


### sofaFit2Grid.m
Conversão de posições fonte para determinadas posições objetivo. Dois método estão disponíveis, 'move': escolha das posições mais proximas entre o grid de entrada e o grid objetivo. 'interp': faz a interpolação de hrirs em posições conhecidas no grid para a projeção em posições objetivo (uso da função interpolateHRTF do matlab audioToolbox).


### sofaResample
Raz o resample de hrirs SOFA para dada taxa de amostragem objetivo. 


###  spec_dist.m 
Calculo da distorção espectral logaritmica entre duas hrirs.


### sph2hor.m 
Conversão entre coordenadas esféricas para horizontais polares


### sph2nav.m
Conversão entre coordenadas esféricas para navegacionais.
