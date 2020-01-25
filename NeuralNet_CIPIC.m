clear all; close all; clc;
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
format short
%% Setting up Feed Forward Neural Network with BP
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']);

load('DADOS_TREINAMENTO\input_target.mat')
[no_PC, no_subjects, no_directions, no_channels] = size(weights);

%% Pr�-processamento INPUT TARGET (normalize)
% [xinpt, sig, mu] = nn_input_preprocess(InputMatrix);

%% Inicializando NET 
[I, ~] = size(InputMatrix_CIPIC);
[O, Ntrn] = size(weights);

H = 2*I+1;                   % No of hidden layers
Nw = (I+1)*H+(H*1)*O;    % No of unknown weights
Ntrneq = Ntrn*O;         % No of training equations
Ndof   = Ntrneq - Nw;    % No of degrees of freedom

%%
net = fitnet([11], 'traincgb'); 
net.trainParam.max_fail = 6;
net.trainParam.showWindow = false;
net.performFcn='msereg';

net.divideFcn= 'dividerand'; % divide the data randomly 
net.divideParam.trainRatio= .7; % we use 70% of the data for training 
net.divideParam.valRatio= .25; % 30% is for validation
net.divideParam.testRatio= .05; % 0% for testing

tol_max = 0.04; % tolerancia m�xima do erro do dataset de testes

%% Treinamento 
disp('Iniciado o treinamento da rede neural')
tic
for channel = 1:no_channels % cada orelha
    x = ((InputMatrix_CIPIC(:, :, channel)));  
    clc;
    for i = 1:no_directions
        disp(i)
        t{i} = ([weights(:, :, i, channel)]);% output target varia com a dire��o 
%         t{i} = mapminmax(weights(:, :, i, channel));%output target varia com a dire��o 
        
        net = configure(net, x, t{i}); %define input e output da rede
        net = init(net);           
        [net_CIPIC{i, channel}, tr_CIPIC{i, channel}] = train(net, x, t{i});
 
 
%       Repete o treinamento para mesma dire��o em caso de perf muito ruim
        cont1 = 0;
         while tr_CIPIC{i, channel}.best_tperf > tol_max
            net = init(net);           
            [net_CIPIC{i, channel}, tr_CIPIC{i, channel}] = train(net, x, t{i});
            cont1 = cont1 + 1;
            if cont1 > 10 % caso n�o consiga atingir a meta em 5 tentativas, utilizar� o ultimo treinamento realizado 
                tr_CIPIC{i, channel}.best_tperf
                break
            end                   
         end   
    end
end
toc

disp('Treinamento completo.')
disp('...')


%% SAVE DATA
save('DADOS_TREINAMENTO\net_treinada_CIPIC.mat','net_CIPIC', 'tr_CIPIC');,...
%                                                 'sig', 'mu');
disp('Dados salvos!')

%% PLOT - MSE for the training set
for k = 1:no_channels
    for l = 1:no_directions
        p(l, k) = tr_CIPIC{l, k}.best_tperf;
    end
end
   
plot(p(:,1)); hold on
plot(p(:,2)); hold off
xlim([0 1250])

legend('L', 'R');
xlabel('Indice da posi��o')
ylabel('RMSE')
title('Performance no mini-grupo de testes')