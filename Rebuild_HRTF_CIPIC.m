clear all; close all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; JUNHO/2019
% Simulação de HRTF individualizada, a partir da RNA treinada e uso de dados
% antropométricos
%% Carregar dados
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']); 
load('DADOS_TREINAMENTO/net_treinada_CIPIC.mat')
load('DADOS_TREINAMENTO/input_target.mat')

%% Definir INPUTs NN
% individuos extraidos: [1:14, 42]
% indivíduos para VALIDAÇÃO: 1('003') e 36('148').
 subject = 1;
% nome do individuo dentro do banco de dados
if subject == 1
    subj = '003';
elseif subject == 36
    subj = '148';
else
    subj = num2str(subject);
end

% separando os 8 parametros necessários para a ann
anthro = importdata('anthro.mat');
% d1, d3 : d6
d1 = anthro.D(subject, [1, 3:6]);      % L
d2 = anthro.D(subject, [9, 11:14]);    % R
d = cat(3, d1, d2);

% x1, x3, x12
x = anthro.X(subject, [1, 3, 12]);

% unir as duas matrizes
InptMtx(:,:,1) = [d(:,:,1), x]';  % L 
InptMtx(:,:,2) = [d(:,:,2), x]';  % R

[no_samples, no_PC, no_directions, no_channels] = size(PC_mtx);

%% Preprocessamento input
% O desvio e a média utilizados são oriundos do dataset de treinamento
% [xinpt, ~, ~] = nn_input_preprocess(InptMtx, sig, mu);

%% Simulação do modelo
% result = cell(no_directions, no_channels);
% 
% disp('Simulação iniciada')
% disp('...')
% % computar valores de saida para dada entrada na ann
% tic
% for n = 1:no_channels
%     for i = 1:no_directions
%         result{i, n} = sim(net_CIPIC{i, n}, (InptMtx(:,:, n)));
%     end
% end
% 
% clc; toc
% disp('Simulação de redes concluida!')

%% SAVE//LOAD DATA
% save(['DADOS_TREINAMENTO\subject', num2str(subject), '_sim.mat'],'result');
% disp('Dados salvos!')

load(['DADOS_TREINAMENTO\subject', num2str(subject), '_sim.mat'])
disp('Dados Carregados')
disp('...')
%% Regeneração da HRTF: obtenção da DTF e HRTF_dir_log ao mesmo tempo
% obtenção da PTF_sim e em sequencia a HRTF_dir_log_sim 
hrtf_dir_log_sim = zeros(no_samples, no_directions, no_channels);

for n = 1:no_channels
    for i = 1:no_directions
        hrtf_dir_log_sim(:, i, n) = PC_mtx(:, :, i, n)*(result{i, n});                                   
        hrtf_dir_log_sim(:, i, n) = hrtf_dir_log_sim(:, i, n) + med_vec2(:, i, n); 
    end
end

%%%% Voltando para HRTF_log10, somando a media das medias extraidas
% [~,no_subjects,~]=size(med_vec);
% %Para estimar a hrtf_log descomente as linhas abaixo:
% hrtf_log_sim = zeros(200, no_directions, no_channels);
% for n = 1:no_channels
%     media = sum(med_vec(:,:,n), 2)/no_subjects;
% %      media = med_vec_ori(:,1,n);
%     for i = 1:no_directions
%         hrtf_log_sim(:, i, n) = hrtf_dir_log_sim(:, i, n) + media;
%     end
% end
% Para calcular a DTF a manipulação acima é desnecessária
hrtf_log_sim = hrtf_dir_log_sim;

% Reshappe para formato CIPIC > [no_azi x no_ele x no_samples]
hrtf_sim = zeros(25, 50, no_samples, no_channels);
dim1 = 25; dim2 = 50;
for i = 1:no_channels
    cont = 0;
    for k = 1:dim1
        for l = 1:dim2
            cont = cont+1;
            hrtf_sim(k,l, :, i) = hrtf_log_sim(:,cont, i);
        end
    end
end

% reshape ITD
% ITD_norm = reshape(ITD_cipic(:,:)', 25, 50, []);

%% CÁlCULO DO ITD
%%% METODO 3: modelo utilizado no dataset CIPIC original, (Woodsworth MODIFIED) 
[new_itd, naz, nel] = itd_metodo3(x(1), x(2));
hrir_final.ITD = new_itd;

%% CORREÇÕES DE FASE: (fase mínima + fase relativa) pelo ITA
clc
for k = 1:dim1 %azimute
    for l = 1:dim2  %elevação 
        hL = 10.^(squeeze(hrtf_sim(k, l, :, 1))./20);
        hR = 10.^(squeeze(hrtf_sim(k, l, :, 2))./20);
        % Aplicando fase minima e de excesso 
        itd = new_itd(k, l);
        [IR_minL, IR_minR] = phase_job(hL, hR, itd, naz(k)); 
        % CHECK fase minima
%          isminphase(IR_minR)
        %save no formato CIPIC
        hrir_final.hrir_l(k, l, :) = IR_minL;
        hrir_final.hrir_r(k, l, :) = IR_minR;
    end
end
hrir_final.name = ['individuo_', subj];


%% Preprocessamento do database para comparação
%Carregando todo o dataset original e obtendo hrtf_Log para todos individuos 
clear hrir_med
onde = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\standard_hrir_database';
path = dir([onde '\subject_*\hrir_final.mat']);
for i = 1 : length( path )
        temph(i) = importdata( fullfile( path(i).folder, path(i).name) );
end
N = length(temph);
for i = 1:N
    hrir(i).L = temph(i).hrir_l; 
    hrir(i).R = temph(i).hrir_r;   
    ITD_medido{i} = temph(i).ITD; 
end
no_samples = length(hrir(1).L);
for i = 1: N
    for m = 1:25
        for n = 1:50
            hrtf(i).logL(m,n,:) = 20*log10(abs(fft(hrir(i).L(m,n,:), no_samples))); 
            hrtf(i).logR(m,n,:) = 20*log10(abs(fft(hrir(i).R(m,n,:), no_samples)));
        end
    end
end
no_subjects = N;
med_vec_ori = zeros(no_samples, no_subjects, no_channels);
for i = 1:no_subjects
   temp1 = zeros(no_samples, 1); 
   temp2 = temp1;
   for j = 1:25
       for k = 1:50
            temp1 = temp1 + squeeze(hrtf(i).logL(j, k, :)) + squeeze(hrtf(i).logR(j, k, :));
       end
   end
   temp1 = temp1/(2*no_directions);
   
   med_vec_ori(:, i) = temp1; % (200 x 30 x 2)
end
for i = 1: N
    for j = 1:25
        for k = 1:50
            hrtf(i).dir_logL(j, k, :) = squeeze(hrtf(i).logL(j, k, :)) - med_vec_ori(:, i);
            hrtf(i).dir_logR(j, k, :) = squeeze(hrtf(i).logR(j, k, :)) - med_vec_ori(:, i);
        end
    end
end
for k = 1:dim1 %azimute
    for l = 1:dim2  %elevação 
        hmedL = 10.^(squeeze(hrtf(subject).dir_logL(k, l, 1:no_samples/2))./20);
        hmedR = 10.^(squeeze(hrtf(subject).dir_logR(k, l, 1:no_samples/2))./20);
        itd2 = temph(subject).ITD(k,l);  
        [hrir_med.hrir_l(k,l,:), hrir_med.hrir_r(k,l,:)] = phase_job(hmedL, hmedR, itd2, naz(k));
    end
end
hrir_med.ITD = temph(subject).ITD;
hrir_med.name = '';
%% Export para SOFA
% SIMULATED DATA
Obj_sim = SOFAconvertCIPIC2SOFA(hrir_final);
Obj_sim.GLOBAL_DatabaseName = ['Individuo', subj];
Obj_sim.GLOBAL_ApplicationName = '';
Obj_sim.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

% MEASURED DATA
Obj_med = SOFAconvertCIPIC2SOFA(hrir_med);
Obj_med.GLOBAL_DatabaseName = ['Individuo', subj];
Obj_med.GLOBAL_ApplicationName = '';
Obj_med.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

save('rebuilt_data_CIPIC_ARI_ITA.mat', 'Obj_sim', 'Obj_med')
disp('all done!')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT - HRTF_log para determinada direção º
clc
% % Azimutes
% naz = [-80, -65, -55 -45:5:45, 55, 65, 80]; 
% % Elevações
% nel = -45 + 5.625*(0:49); 
fs = 44100;
Nt = no_samples;  
% f_sw = linspace(0,(Nt-1)*fs/Nt, Nt);
freq = linspace(0, fs-fs/Nt, Nt);

% DEFINA A DIREÇÃO DA HRTF
%  1 6 13 20 25
azz = 25; % azimute 1:25 % 13 representa 0º de azimute
ell = 9; % elevação 1:50 % 9 representa 0º de elevação e 41 - 180°

%%% Medido %%%
% L
gL = squeeze(hrir_med.hrir_l(azz,ell,:));
g_smoothL = 20*log10(abs(fft(gL)));
% g_smoothL = smooth(gLf, 'lowess');
% R 
gR = squeeze(hrir_med.hrir_r(azz,ell,:));
g_smoothR = 20*log10(abs(fft(gR)));
% g_smoothR = smooth(gRf, 'lowess');

%%% Simulado %%%
% hori = (squeeze(hrtf_sim(azz, ell, 1:Nt/2, 1)));
% h_prop_smoothL = smooth(hori, 'lowess');

% Apos de fase
% L
hrtf_propL = squeeze(hrir_final.hrir_l(azz, ell, :));
h_prop_smoothL = 20*log10(abs(fft(hrtf_propL)));
% h_prop_smoothL = (smooth(h_propL, 'lowess'));
% R 
hrtf_propR = squeeze(hrir_final.hrir_r(azz, ell, :));
h_prop_smoothR = 20*log10(abs(fft(hrtf_propR)));
% h_prop_smoothR = (smooth(h_propR, 'lowess'));

%%%% PLOT %%%%
figure()
semilogx(freq(1:Nt/2), g_smoothL(1:Nt/2), 'lineWidth', 1.5,  'Color', 'blue'); hold on
semilogx(freq(1:Nt/2), g_smoothR(1:Nt/2), 'lineWidth', 1.5, 'Color', 'red'); 


semilogx(freq(1:Nt/2), h_prop_smoothL(1:Nt/2), '-.r', 'lineWidth', 1.5, 'Color', 'blue'); 
semilogx(freq(1:Nt/2), h_prop_smoothR(1:Nt/2), '-.r', 'lineWidth', 1.5, 'Color', 'red'); hold off

legend('Original Esquerda', 'Original Direita', ... 
       'Proposta Esquerda', 'Proposta Direita', ...
               'Location', 'best', 'lineWidth', 1.5)
    
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');

azim = num2str(naz(azz));
elee = num2str(nel(ell));
title(['DTF Individuo ', subj,', azimute ' azim 'º, elevação ' elee 'º.' ]);

grid on 
% ylim([-20 20])
xlim([200, 2e4])
set(gca,'FontSize',13)


%% PLOT - HRIR [Simulada vs Medida]
figure()
azi = 1;
ele = 9;

%%% ESQUERDA %%%
subplot(2,1,1)
% Medida
plot( real(squeeze(hrir_med.hrir_l(azi,ele,:))), 'linewidth', 1.5);hold on
% Proposta
plot(real(squeeze(hrir_final.hrir_l(azi, ele, :))), '-.r', 'linewidth', 1.3); hold off

xlabel('Amostras')
ylabel('Esquerda')
% ylim([-1 1])
title(['HRIR Invididuo ', num2str(subj), ', azimute ', num2str(naz(azi)),...
                                   '°, elevação ', num2str(nel(ele)), '°']);
legend('Original', 'Simulada')
set(gca,'FontSize',12)
axis tight

%%% DIREITA %%%
subplot(2,1,2)
plot( real(squeeze(hrir_med.hrir_r(azi,ele,:))), 'linewidth', 1.5); hold on
plot( real(squeeze(hrir_final.hrir_r(azi, ele, :))) , '-.r','linewidth', 1.3); hold off

xlabel('Amostras')
ylabel('Direita')
legend('Original', 'Simulada')

% ylim([-1 1])
set(gca,'FontSize',12)
axis tight






