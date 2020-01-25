% clear all; clc
% % DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Novembro/2019
% %% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas')
% load('Resultados/deepnet_treinada.mat')
% addpath([pwd, '/Functions'])

%% Dados ANTHROPOMETRICOS (INPUT DATA)
% Carregar dados 
anthro = importdata('anthro.mat');
subj = 1;
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

%% Normalização do input
inpt_subject = nn_input_preprocess(anthro_data(:, subj,:));
input = inpt_subject(:,subj,:);
% input = anthro_data(:, subj,:);

%% Simulação do modelo
% [dim1, dim2, no_channels] = size(trainedNet);
% no_samples = 200;
% result = zeros(dim1, dim2, no_samples, no_channels);


% disp('Iniciando simulação de Redes ...')
% tic
% wait = waitbar(0,'Processando HRTFs...');
% for k = 1:dim1
%     for l = 1:dim2
% %         for m = 1:no_channels
%             result(k, l, :, m) = ((predict(trainedNet{k,l,m}, input(:,:,m)))*4)-2;
% %         end
%     end
%     waitbar(k/dim1, wait)
% end
% toc
% close(wait)


net = importKerasNetwork('Python/keras_fitted.h5');
input = permute(input, [2, 1 , 3]);
result = ((predict(net, input(1,:,1)))*4)-2;
disp('Concluido!')


%% Preprocessamento do database para comparação
%Carregando todo o dataset original e obtendo hrtf_Log para todos individuos 

path = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_*\hrir_final.mat');
for jj = 1 : subj
        data(jj) = importdata( fullfile( path(jj).folder, path(jj).name) );
end
N = length(data);
for jj = 1:N
    hrir(jj).L = data(jj).hrir_l; 
    hrir(jj).R = data(jj).hrir_r;   
end
hrir_final_measured = data(subj);


%% Export SOFA
% hrir_final.hrir_l = result(:,:,:,1);
% hrir_final.hrir_r = result(:,:,:,2);
% hrir_final.name = ['subject', num2str(subj)];
% % simulado
% Obj_sim = SOFAconvertCIPIC2SOFA(hrir_final);
% Obj_sim.GLOBAL_DatabaseName = ['Individuo', subj];
% Obj_sim.GLOBAL_ApplicationName = '';
% Obj_sim.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
% % medido
% Obj_med = SOFAconvertCIPIC2SOFA(hrir_final_measured);
% Obj_med.GLOBAL_ApplicationName = '';
% Obj_med.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
% 
% %%%% export SOFA file %%%%
% file_path = fullfile(['HRIRs_sim\deepnet_individuo_' num2str(subj) '.sofa']);
% disp(['Saving:  ' file_path])
% compression = 0;
% SOFAsave(file_path, Obj_sim, compression);
% 
% 
% %%%% save to analysis %%%%
% save('rebuilt_deep_CIPIC.mat', 'Obj_sim', 'Obj_med')
% disp('all done!')
% toc


%% Plot 

azi = 13; % AZIMUTE
ele = 9; % ELEVACAO

y0 = squeeze(hrir(subj).L(azi, ele, :));
% y1 = squeeze(result(azi, ele, :, 1));
y1 = result;
no_samples = 200;
% y0 = zscore(squeeze(hrir(1).L(azi, ele, :)),1,'all'); % medido
% y1 = zscore(squeeze(result(azi, ele, :, 1)),1,'all'); % simulado

% TEMPO 
figure()
plot(y0,'linewidth', 2); hold on
plot(y1, '-.r','linewidth', 1.5)
legend('Medido', 'Simulado', 'Location', 'best')
xlabel('Amostras')
ylabel('Amplitude')

% FREQ
fs = 44100; 
N = no_samples;
freq = linspace(0, fs-fs/N, N);


H0 = 20*log10(abs(fft(y0)));
H1 = 20*log10(abs(fft(y1)));

figure()
semilogx(freq(1:no_samples/2), H0(1:no_samples/2), 'linewidth', 2); hold on 
semilogx(freq(1:no_samples/2), H1(1:no_samples/2), '-.r','linewidth', 1.5);
legend('Medido', 'Simulado', 'Location', 'best')
title('Magnitude')
xlabel('Frequência [Hz]')
ylabel('Amplitude [dB]')
grid on 
axis tight



