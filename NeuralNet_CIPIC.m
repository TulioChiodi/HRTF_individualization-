clear all; close all; clc;
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019

%% Setting up Feed Forward Neural Network with BP
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'

load('DADOS_TREINAMENTO\input_target.mat')
[no_PC, no_individuos, no_directions, no_channels] = size(weights);
save('DADOS_TREINAMENTO\python_convert.mat', 'InputMatrix', 'ITD_cipic', 'weights');
net = fitnet(16); %cria��o da rede com 12 n�s no hidden layer
net.trainFcn = 'trainbr'; % fun��o de treinamento backwardpropagation
net.trainParam.goal = 1e-6;   % error goal
net.trainParam.epochs = 2000; % numero de itera��es
tol_max = 0.09; % tolerancia m�xima do erro do dataset de testes

net.trainParam.max_fail = 20;

% bayesopt2
% lstm com hrir
%% Treinamento 
disp('Iniciado o treinamento da rede neural')

tic

% Uma network para cada dire��o// (� a forma mais eficiente??)
for channel = 1:no_channels % cada orelha
    x = InputMatrix(:, :, channel);   
    clc;
    for i = 1:no_directions

        t{i} = [weights(:, :, i, channel)];% output target varia com a dire��o 

%         t{i} = [weights(:, :, i, channel); ITD_cipic(i, :)];% output target varia com a dire��o 
        net = configure(net, x, t{i}); %define input e output da rede
        net = init(net);
        [net_CIPIC{i, channel}, tr_CIPIC{i, channel}] = train(net, x, t{i});
        
        %repete o treinamento para mesma dire��o em caso de perf muito ruim
        cont1 = 0;
         while tr_CIPIC{i, channel}.best_tperf > tol_max
            net = init(net);           
            [net_CIPIC{i, channel}, tr_CIPIC{i, channel}] = train(net, x, t{i});
%             disp(['repeti��o para itera��o: ', num2str(i)])
            cont1 = cont1 + 1;
            if cont1 > 5 % caso n�o consiga atingir a meta em 5 tentativas, utilizar� o ultimo treinamento realizado 
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
save('DADOS_TREINAMENTO\net_treinada_CIPIC.mat','net_CIPIC', 'tr_CIPIC');
disp('Dados salvos!')

%% avalia�� de performance

for k = 1:no_channels
    for l = 1:no_directions
        p(l, k) = tr_CIPIC{l, k}.best_tperf;
    end
end
   
plot(p(:,1)); hold on
plot(p(:,2)); hold off
ylim([0 3]);
xlim([0 1250])
legend('L', 'R');