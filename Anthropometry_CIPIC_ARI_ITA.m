clear all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Fevereiro/2020
%% GENERAL INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~ ASSEMBLE ANTHROPOMETRY FROM CIPC, ARI AND ITS DATABASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']);
% CIPIC 
anthro_CIPIC = load('anthro_CIPIC.mat');

% ARI 
anthro_ARI = load('anthro_ARI.mat');

% ITA
anthro_ITA = load('anthro_ITA.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CIPIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orelhas
d1_C = anthro_CIPIC.D(:, 1:8);
d2_C = anthro_CIPIC.D(:, 9:16);
% d3_C = anthro_CIPIC.D;
% t1_C = anthro_CIPIC.theta(:, 1:2);
% t2_C = anthro_CIPIC.theta(:, 3:4);
% t3_C = anthro_CIPIC.theta;

% Cabeca e torso
x_C = anthro_CIPIC.X(:,[1,3]);

% Concatena
data_CIPIC(:,:,1) = [x_C, d1_C].'; 
data_CIPIC(:,:,2) = [x_C, d2_C].'; 
% data_CIPIC(:,:,1) = [x_C, t3_C, d3_C].';
% data_CIPIC(:,:,2) = [x_C, t3_C, d3_C].';

%%% Remorção de indivíduos sem antropometria
for k = 1:2
    data_in = data_CIPIC(:,:,k);
    remove_CIPIC = any(isnan(data_in),1); % identifica quais colunas não possuem pelo menos uma das medidas necessárias 
    remove_CIPIC(1) = 1; % remover individuo pra posterior teste
    data_in(:,remove_CIPIC) = [];
    CIPIC_anthro(:,:,k) = data_in;    
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orelhas
d1_A = anthro_ARI.D(:, 1:8);
d2_A = anthro_ARI.D(:, 9:16);
% d3_A = anthro_ARI.D(:,1:16);
% t1_A = anthro_ARI.theta(:, 1:2);
% t2_A = anthro_ARI.theta(:, 3:4);
% t3_A = anthro_ARI.theta;

% Cabeca e torso
x_A = anthro_ARI.X(:,[1,3]);

% Concatena
data_ARI(:,:,1) = [x_A, d1_A].'; 
data_ARI(:,:,2) = [x_A, d2_A].'; 
% data_ARI(:,:,1) = [x_A, t3_A, d3_A].'; 
% data_ARI(:,:,2) = [x_A, t3_A, d3_A].'; 

%%% Remorção de indivíduos sem dados anthropométricos %%%
for m = 1:2
    data_in = data_ARI(:,:,m);
    remove_ARI = any(isnan(data_in),1); % identifica quais colunas não possuem pelo menos uma das medidas necessárias 
    data_in(:,remove_ARI) = [];
    ARI_anthro(:,:,m) = data_in;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_I = anthro_ITA.D;
x_I = anthro_ITA.X;
ITA_anthro = [x_I; d_I]/10;
ITA_anthro(:,:,2) = ITA_anthro;


%% EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anthro = cat(2, [CIPIC_anthro, ARI_anthro, ITA_anthro]);

% imagesc(anthro(:,:,2))

%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save('DADOS_TREINAMENTO/input_CIPIC.mat', 'CIPIC_anthro')
% save('DADOS_TREINAMENTO/input_CIPIC_ARI.mat', 'anthro')
save('DADOS_TREINAMENTO/input_CIPIC_ARI_ITA.mat', 'anthro')

% save('remove_ARI.mat', 'remove_ARI')
% save('remove_CIPIC.mat', 'remove_CIPIC')
% disp('Dados Salvos!')