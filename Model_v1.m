clear all; close all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; JUNHO/2019
%%
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']);
load('DADOS_TREINAMENTO/net_treinada_CIPIC.mat')
load('DADOS_TREINAMENTO/input_target.mat')
tic
%% Defnindo INPUT
% individuos extraidos: [1:14, 42]
% indivíduos para VALIDAÇÃO: 1('003') e 36('148').

% subject = input('Escolha um numero para ID: ');
subject = 140;
% nome do individuo dentro do banco de dados
if subject == 1 || subject == 36 
%      error('Esta ID está em uso!')
subj = '003';
else
    subj = num2str(subject);
end

% separando os 8 parametros necessários para a ann
% d1, d3 : d6
% dL = anthro.D(subject, [1, 3:6]);      % L
% dR = anthro.D(subject, [9, 11:14]);    % R
dL = [1.45, 1.6, 1.6, 6.54, 3.83];
dR = [1.6, 1.6, 1.7, 6.4, 3.5];
d = cat(3, dL, dR);

% x1, x3, x12
% x1 = head width
% x3 = head depth
% x12 = shoulder width
% x = anthro.X(subject, [1, 3, 12]);
x = [20, 27, 49];

% unir as duas matrizes
InptMtx(:,:,1) = abs([d(:,:,1), x]');  % L 
InptMtx(:,:,2) = abs([d(:,:,2), x]');  % R


%% Simulação do modelo
[no_samples, no_PC, no_directions, no_channels] = size(PC_mtx);
result = cell(no_directions, no_channels);

disp('Simulação iniciada')
disp('...')
% computar valores de saida para dada entrada na ann
tic
for n = 1:no_channels
    for i = 1:no_directions
        result{i, n} = sim(net_CIPIC{i, n}, (InptMtx(:,:, n)));
    end
end

clc; toc
disp('Simulação de redes concluida!')


%% SAVE//LOAD DATA
% save(['DADOS_TREINAMENTO\subject', num2str(subject), '_sim.mat'],'result', 'ITD');
% disp('Dados salvos!')
% 
% load(['DADOS_TREINAMENTO\subject', num2str(subject), '_sim.mat'])
% disp('Dados Carregados')

disp('...')
%% Regeneração da DTF (PTF)
% obtenção da PTF_sim e em sequencia a HRTF_dir_log_sim 
DTF_sim = zeros(no_samples, no_directions, no_channels);

for n = 1:no_channels
    for i = 1:no_directions
        DTF_sim(:, i, n) = PC_mtx(:, :, i, n)*(result{i, n});                                   
        DTF_sim(:, i, n) = DTF_sim(:, i, n) + med_vec2(:, i, n); 
    end
end


%% Voltando para HRTF_log10, somando a media das medias extraidas
% Reshappe para formato CIPIC > [no_azi x no_ele x no_samples]
dim1 = 25; dim2 = 50;
hrtf_sim = zeros(dim1, dim2, no_samples, no_channels);
for i = 1:no_channels
    cont = 0;
    for k = 1:dim1
        for l = 1:dim2
            cont = cont+1;
            hrtf_sim(k,l, :, i) = DTF_sim(:,cont, i);
        end
    end
end

%% CÁlCULO DO ITD
%%% Modelo utilizado no dataset CIPIC original, (Woodsworth MODIFIED) 
[new_itd, naz, nel] = itd_metodo3(x(1), x(2));
hrir_final.ITD = new_itd;

%% CORREÇÕES DE FASE: (fase mínima + fase relativa) pelo ITA
for k = 1:dim1 %azimute
    for l = 1:dim2  %elevação 
        hL = 10.^(squeeze(hrtf_sim(k, l, :, 1))./20);
        hR = 10.^(squeeze(hrtf_sim(k, l, :, 2))./20);
        % Aplicando fase minima e de excesso 
        itd = new_itd(k, l);
        [IR_minL, IR_minR] = phase_job(hL, hR, itd, naz(k)); 
        %save no formato CIPIC
        hrir_final.hrir_l(k, l, :) = IR_minL;
        hrir_final.hrir_r(k, l, :) = IR_minR;
    end
end
hrir_final.name = ['individuo_', subj];


%% Export para SOFA
Obj = SOFAconvertCIPIC2SOFA(hrir_final);
Obj.GLOBAL_DatabaseName = ['Individuo', subj];
Obj.GLOBAL_ApplicationName = '';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');

%%%% save SOFA file %%%%
file_path = fullfile(['HRIRs_sim\individuo_' subj '.sofa']);
disp(['Saving:  ' file_path])
compression = 0;
SOFAsave(file_path, Obj, compression);

disp('all done!')
toc

%% DAFF export 
% hrtf = itaHRTF('sofa', file_path);
% hrtf.writeDAFFFile('HRIRs_sim\')

%% SAVE CIPIC DATA
% save(['HRIRs_sim\individuo_' subj '.mat'], 'hrir_final');
% toc
% disp('Dados salvos!')

%% Plot geometry mesh
% Obj.SourceView_Type = 'spherical';
% Obj.API.Dimensions.SourceView  = 'MC';
% Obj.SourceView = Obj.SourcePosition;
% SOFAplotGeometry(Obj)

%% PLOT - HRIR
figure()

azz = 25; % azimute 1:25 % 13 representa 0º de azimute
ell = 9; % elevação 1:50 % 9 representa 0º de elevação

subplot(2,1,1)
plot( real(squeeze(hrir_final.hrir_l(azz, ell, :))), 'linewidth', 1.3);
xlabel('Amostras')
ylabel('Esquerda')
title(['HRIR Invididuo ', num2str(subj), ', azimute ', num2str(naz(azz)),...
                                   '°, elevação ', num2str(nel(ell)), '°']);
set(gca,'FontSize',12)
axis tight

subplot(2,1,2)
plot( real(squeeze(hrir_final.hrir_r(azz, ell, :))) , 'linewidth', 1.3);
xlabel('Amostras')
ylabel('Direita')
set(gca,'FontSize',12)
axis tight


%% PLOT - DTF / HRTF
fs = 44100;
Nt = no_samples; 
f_sw = linspace(0, fs-fs/Nt, Nt);
K = 1:Nt/2; 

azz = 13;
ell = 9; 


figure()

%Smoothie  
smoo_func = 'sgolay';

% Simulado
h = 20*log10(abs(fft(squeeze(hrir_final.hrir_l(azz, ell, :)))));
% h = h - mean(h);
h_smooth = smooth(h, smoo_func);

g = 20*log10(abs(fft(squeeze(hrir_final.hrir_r(azz, ell, :)))));
% g = g - mean(g);
g_smooth = smooth(g, smoo_func);

% PLOT
semilogx(f_sw(K), h_smooth(K), 'lineWidth', 1.6, 'Color', 'red'); hold on
semilogx(f_sw(K), g_smooth(K), 'lineWidth', 1.6, 'Color', 'blue'); hold off

legend('Esquerda', 'Direita', 'Location', 'best',...
                                        'lineWidth', 1.5)
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
azi = num2str(naz(azz));
ele = num2str(nel(ell));
title([ 'azimute ' num2str(naz(azz)) 'º, elevação ' num2str(nel(ell)) 'º.' ]);
grid on 
% ylim([-20 20])
xlim([200, 2e4])
set(gca,'FontSize',13)
