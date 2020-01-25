clear all; clc 
%%% INDIVIDUALIZED HRTF USING ANTHROPOMETRIC DATA AND DEEP NEURAL NETWORKS
% Davi Carvalho/ novembro 2019
addpath([pwd, '/Functions'])

%% Dados ANTHROPOMETRICOS (INPUT DATA)
% Carregar dados 
anthro = importdata('anthro.mat');

% separando os 8 parametros necessários para a ann
% ear dimensions
d1 = anthro.D(:, 1:8);
d2 = anthro.D(:, 9:16);
d3 = anthro.D(:,1:16);

t1 = anthro.theta(:, 1:2);
t2 = anthro.theta(:, 3:4);
t3 = anthro.theta(:, 1:4);
% head and torso
x = anthro.X;
% concatena
anthro_data(:,:,1) = [x, t3, d3].';
anthro_data(:,:,2) = [x, t3, d3].';

% anthro_data(:,:,1) = [x, t1, d1].'; 
% anthro_data(:,:,2) = [x, t2, d2].'; 
% 
%% Remorção de indivíduos sem dados anthropométricos
for k = 1:2
    data_in = anthro_data(:,:,k);
    remove = any(isnan(data_in),1); % identifica quais colunas não possuem pelo menos uma das medidas necessárias 
    remove(1) = 1; % remover individuo pra posterior teste
    data_in(:,remove) = [];
    anthro_dataOK(:,:,k) = data_in;    
end 

%% Normalização do input
[input, sig, mu] = nn_input_preprocess(anthro_dataOK);
% surf(input(:,:,2))
% view([90 0])

%% LOAD  TARGET DATA: Carregar o banco de HRIRs  
% Defina como path o local onde estão as pastas com as pastas dos arquivos .mat do banco de dados.
diretorio = (['B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\standard_hrir_database']);
path = dir([diretorio '\subject_*\hrir_final.mat']);
% Importa todo o banco de indivíduos para 'data'
for jj = 1 : length( path )
    data(jj) = importdata( fullfile( path(jj).folder, path(jj).name) );
end

% Remove indivíduos sem dados antropométricos e indivíduos teste ann
fields = fieldnames(data); %salva nome dos fields do struct para reconstrução 
data  = struct2cell(data); %transforma em celula para remover indivíduos
data(:, :, remove) = []; %remove indivíduos 
data = cell2struct(data, fields, 1); %volta para struct

%% Extração de dados [HRIR]
N = length(data);
for jj = 1:N
    hrir(jj).L = data(jj).hrir_l; 
    hrir(jj).R = data(jj).hrir_r;
%     target(:,:,:, jj, 1) = 20*log10(abs(fft(hrir(jj).L, 200, 3)))/32;
%     target(:,:,:, jj, 2) = 20*log10(abs(fft(hrir(jj).R, 200, 3)))/32;
    target(:,:,:, jj, 1) = (hrir(jj).L + 2)./4;
    target(:,:,:, jj, 2) = (hrir(jj).R + 2)./4;
end

%python transfer
save('Python/data_for_train.mat', 'target', 'input', 'sig', 'mu') 

%% TREINAMENTO
[no_azi, no_ele, no_samples,no_subjects, no_channels] =  size(target);
trainedNet = cell(no_azi, no_ele, no_channels);
probability = 0.7;
nodes = 64;
tic
disp('Iniciado Treinamento,')
for k = 1:no_channels
%%% ARQUITETURA DA REDE %%%
    % Camadas
    layers = [
        imageInputLayer([1,37,1])
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot')
        reluLayer
        dropoutLayer(probability)        
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot')
        reluLayer
        dropoutLayer(probability)         
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot')
        reluLayer
        dropoutLayer(probability)
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot')
        reluLayer
        dropoutLayer(probability)
        fullyConnectedLayer(nodes,'WeightsInitializer','glorot')
        reluLayer
        dropoutLayer(probability)
        fullyConnectedLayer(200)
        reluLayer
        regressionLayer
        ];
                    
    for l = 1:no_azi
        disp(['Azimute: ', num2str(l),' de 25 '] )
        tic       
        for m = 1:no_ele
            disp(['Elevação ', num2str(m),' de 50 ...'] )        
            
            % SPLIT DATASETS for validation
            X = input(:,:,k);
            Y = squeeze(target(l,m,:,:,k));
             trainRatio = 0.9;
             valRatio = 0.1;
             testRatio = 0;
            [XTrain, XVal, ~] = divideint(X,trainRatio,valRatio,testRatio);          
            [YTrain, YVal, ~] = divideint(Y,trainRatio,valRatio,testRatio);          
            XTraiin(1,:,1,:)  = XTrain;
            YTraiin(1,1,:,:)  = YTrain;
            XValidation(1,:,1,:)  = XVal;
            YValidation(1,1,:,:)  = YVal;
            
            %%% Opções de treinamento %%%           
            if m == 1
            options = trainingOptions('adam', ...
                          'MaxEpochs', 20000, ...
                          'InitialLearnRate', 0.001,... 
                          'GradientDecayFactor', 0.9,...
                          'SquaredGradientDecayFactor', 0.999,...
                          'ValidationData', {XValidation, YValidation}, ...
                          'Shuffle', 'once', ...
                          'Verbose', 1);...                              
%                               'Plot','training-progress'); %mostra grafico

            else % Transfer Learning da posição anterior
                layers = trainedNet{l,(m-1),k}.Layers; 
                options = trainingOptions('adam', ...
                          'InitialLearnRate', 0.001,... 
                          'GradientDecayFactor', 0.9,...
                          'SquaredGradientDecayFactor', 0.999,...
                          'MaxEpochs', 1000, ...
                          'ValidationData', {XValidation, YValidation}, ...
                          'Shuffle', 'once', ...
                          'Verbose', 0);  

            end            
            %%% TREINAMENTO %%%
            trainedNet{l,m,k} = trainNetwork(XTraiin, YTraiin, layers, options);
        end
        toc
    end
end
toc          
            
%% SAVE DATA 
save('Resultados/deepnet_treinada.mat', 'trainedNet', 'sig', 'mu');
disp('Dados Salvos!')

%% Desligar pc
% system('shutdown -s')
% disp('O sistema será desligado em 60s!')











