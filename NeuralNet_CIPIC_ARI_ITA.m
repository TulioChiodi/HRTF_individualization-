clear all; close all; clc;
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; SETEMBRO/2019
format short
%% LOAD
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
% inputs
load('DADOS_TREINAMENTO/input_target.mat'); clear weights med_vec  ITD_cipic
load('DADOS_TREINAMENTO/input_matrix_ARI.mat'); clear med_vec2
load('DADOS_TREINAMENTO/input_ITA.mat')
% target
load('DADOS_TREINAMENTO/target_cipic_ari_ita.mat') %weights

%% Unir inputs e targets (ARI+CIPIC)
InputMatrix_CIPIC = InputMatrix_CIPIC(1:7,:,:);
InputMatrix_ARI = InputMatrix_ARI(1:7,:,:);
InputMatrix_ITA = InputMatrix_ITA;
InptMtx = (cat(2, [InputMatrix_CIPIC, InputMatrix_ARI, InputMatrix_ITA]));

%% Pr�-processamento INPUT (normalize)
% [xinpt, sig, mu] = nn_input_preprocess(InptMtx);

%% Setting up Feed Forward Neural Network with BP
net = fitnet([12]); %cria��o da rede com 12 n�s no hidden layer
% fun��o de treinamento
net.trainFcn = 'traincgb';  %cgb
net.trainParam.epochs = 500; % numero de itera��es
tol_max = 0.025; % tolerancia m�xima do erro para teste subset 
net.trainParam.max_fail = 10;

net.divideFcn= 'dividerand'; % divide the data randomly 
net.divideParam.trainRatio= .7; % we use 70% of the data for training 
net.divideParam.valRatio= .25; % 30% is for validation
net.divideParam.testRatio= .05; % 10% for testing

net.trainParam.showWindow = false;

%% Treinamento 
disp('Iniciado o treinamento da rede neural')
[no_PC, no_subjects, no_directions, no_channels] = size(weights);

tic
temp1 = 0;
for channel = 1:no_channels % cada orelha  
    for i = 1:no_directions
        disp(i)
        % definir input para cada canal
        x = (InptMtx(:,:, channel));
        %definir target para cada dire��o, e mapeamento entre 0 e 1
        t{i} = ([weights(:, :, i, channel)]);
         
        % treinamento 
        net = init(net);
        [net_CIPIC_ARI_ITA{i, channel}, tr_CIPIC_ARI_ITA{i, channel}] = train(net, x, t{i});
        
%         %caso de perf muito ruim, atualiza��o de pesos e novo treinamento 
        cont = 0;
             while tr_CIPIC_ARI_ITA{i, channel}.best_tperf > tol_max
                net = init(net);           
                [net_CIPIC_ARI_ITA{i, channel}, tr_CIPIC_ARI_ITA{i, channel}] = train(net, x, t{i});
                cont = cont + 1;
                if cont > 5 % teste de 5 tentativas para atingir a meta
                             tr_CIPIC_ARI_ITA{i, channel}.best_tperf
                    break
                end
             end         
    end
end
toc

disp('Treinamento completo.')
disp('...')

%% SAVE DATA
save('DADOS_TREINAMENTO/net_treinada_CIPIC_ARI_ITA.mat','net_CIPIC_ARI_ITA', 'tr_CIPIC_ARI_ITA');
disp('Dados salvos!')


%% PLOT: performance do grupo de treinamento
for i = 1:no_directions
    
    pl(i) = tr_CIPIC_ARI_ITA{i, 1}.best_tperf;
    pl2(i) = tr_CIPIC_ARI_ITA{i, 2}.best_tperf;
end
plot(pl); hold on 
plot(pl2);


xlim([0 1250])
ylim([0 0.07])

legend('L', 'R');
xlabel('Indice da posi��o')
ylabel('RMSE')
title('Performance no mini-grupo de testes')
