clc; close all 
% % DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Novembro/2019
% Simulação de HRTFS com rede treinada em Python (Keras - tensorflow backend)
% %% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas')
addpath([pwd, '/Functions'])

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
input = permute(input, [2, 1 , 3]);

%% Carregar do modelo
no_azi = 25; no_ele = 50; no_channels = 2; no_samples = 200;
trainedNet = cell(no_azi, no_ele, no_channels);

disp('Carregando Modelo ...')
tic
wait = waitbar(0, 'Carregando Modelo ...');
for k = 1:no_azi
    for l = 1:no_ele
        for m = 1:no_channels
            % Path com as networks treinadas
            path_net = [pwd '/Python/Deep_Keras_C_64/net_', ...
                         num2str(k) '_' num2str(l) '_' num2str(m) '.h5'];
            trainedNet{k,l,m} = importKerasNetwork(path_net);
        end
    end
    waitbar(k/no_azi, wait)
end
toc
close(wait)

%% Simulacao 
result = zeros(no_azi, no_ele, no_samples, no_channels);
disp('Iniciando simulação de Redes ...')
tic
wait = waitbar(0,'Processando HRTFs...');
for k = 1:no_azi
    for l = 1:no_ele
        for m = 1:no_channels
            result(k, l, :, m) = ((predict(trainedNet{k,l,m}, input(:,:,m)))*32);
        end
    end
    waitbar(k/no_azi, wait)
end
toc
close(wait)

disp('Concluido!')


%% Preprocessamento do database para comparação
%Carregando todo o dataset original e obtendo hrtf_Log para todos individuos 

path = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_*\hrir_final.mat');
for jj = 1 : subj
    data(jj) = load( fullfile( path(jj).folder, path(jj).name) );
end
N = length(data);
for jj = 1:N
    hrir(jj).L = data(jj).hrir_l; 
    hrir(jj).R = data(jj).hrir_r;   
end

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
close all 
azi = 1; % AZIMUTE
ele = 9; % ELEVACAO
% 
y0 = squeeze(hrir(subj).L(azi, ele, :));
y1 = squeeze(result(azi, ele, :, 1));
no_samples = 200;
% y0 = zscore(squeeze(hrir(subj).L(azi, ele, :)),1,'all'); % medido
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



