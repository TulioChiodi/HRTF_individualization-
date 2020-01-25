clear all; close all; clc
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; JUNHO/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Avaliacao da HRTF simulada, por meio de comparacoes entre HRTF medida,
% simulada e generica
%% PATHs
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']);

% Cabecas 
addpath('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Cabecas');


%% LOAD HRTFs SIMULADA e MEDIDADA
% load('rebuilt_data_CIPIC.mat')
load('rebuilt_data_CIPIC_ARI_ITA.mat')
dtf_med = fit_2_CIPIC_grid(Obj_med); 
dtf_sim = fit_2_CIPIC_grid(Obj_sim); 

fs = dtf_med.Data.SamplingRate;
no_samples = size(dtf_med.Data.IR,3);

%% LOAD Generic HRTFs
%Fabian
generic_head = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa');

% MIT KEMAr large pinna
% generic_head = SOFAload('mit_kemar_large_pinna.sofa');

% MIT 
% generic_head = SOFAload('mit_kemar_normal_pinna.sofa');


%%% DTF %%% 
generic_fitted = fit_2_CIPIC_grid(generic_head, 'interp');
[dtf_gen, ~] = SOFAhrtf2dtf(generic_fitted, 'log');

%% LSD
% Calculo do erro
for k = 1:length(dtf_med.Data.IR)     
    % separando IR
    hrir_med = squeeze(dtf_med.Data.IR(k, 2, :));
    hrir_sim = squeeze(dtf_sim.Data.IR(k, 2, :));
    hrir_gen = squeeze(dtf_gen.Data.IR(k, 2, :));
    
    %%% LSD %%%
    sd_sim(k) = spec_dist(hrir_sim, hrir_med, fs, 200, 1.3e4);
    sd_gen(k) = spec_dist(hrir_sim, hrir_gen, fs, 200, 1.3e4);
    
    %%% RMSE %%%
    rmse_sim(k) = 20*log10(rms(abs(hrir_med - hrir_sim)));
    rmse_gen(k) = 20*log10(rms(abs(hrir_med - hrir_gen(1:no_samples))));
end

%% Mapa projeção cilindrica
% Simulação 
figure()
pos_sph = round(dtf_med.SourcePosition);
scatter(pos_sph(:,1), pos_sph(:,2), 30, sd_sim, 'filled', 'square')
title('Distorção espectral logaritmica - HRTFs Simuladas')
xlabel('Azimute')
ylabel('Elevação')
axis tight
c = colorbar; caxis([5 12]); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)

% Generica
figure()
pos_sph = round(dtf_med.SourcePosition);
pl_gen = scatter(pos_sph(:,1), pos_sph(:,2), 30, sd_gen, 'filled', 'square');
title('Distorção espectral logaritmica - HRTFs Genéricas')
xlabel('Azimute')
ylabel('Elevação')
axis tight
c = colorbar; caxis([5 12]); colormap jet
c.Label.String = 'Distorção Espectral [dB]';
set(gca,'FontSize',12)


%% Plot  HRTFs SURFACE
figure()
type = 'MagMedian';
subplot(211)
SOFAplotHRTF(Obj_med, type);
title('HRTFs medidas')
axis tight

subplot(212)
SOFAplotHRTF(Obj_sim, type);
title('HRTFs simuladas')
axis tight
colormap jet


%% LSD - Log Spectral Distance
% % ELEVAï¿½ï¿½O ESPECï¿½FICA
% 
% elevation = 0;
% freq_max = 1.3e4;
% 
% % encontrar posiï¿½ï¿½es relativas a elevaï¿½ï¿½o escolhida
% idxmed = find(round(dtf_med.SourcePosition(:,2)) == elevation);
% idxmed = sort(idxmed);
% 
% for k = 1:length(idxmed) 
%     hrir_med = squeeze(dtf_med.Data.IR(idxmed(k),1,:));
%     hrir_sim = squeeze(dtf_sim.Data.IR(idxmed(k),1,:));
%     hrir_gen = squeeze(dtf_gen.Data.IR(idxmed(k),1,:));
%     %%% LSD %%%
%     sd_sim(k) = spec_dist(hrir_sim,hrir_med, fs, 200, freq_max);
%     sd_gen(k) = spec_dist(hrir_sim,hrir_gen, fs, 200, freq_max);    
%     
% %     %%% RMSE %%%
% %     hrir_msrd = normalize(hrir_msrd, 'range');
% %     hrir_sim = normalize(hrir_sim, 'range');
% %     hrir_fabian = normalize(hrir_fabian, 'range');
%     
%     rmse_sim(k) = 20*log10(rms(abs(hrir_med - hrir_sim)));
%     rmse_gen(k) = 20*log10(rms(abs(hrir_med - hrir_gen(1:no_samples))));
% end
% %% PLOT
% %%% LSD %%%
% figure()
% hold on 
% bar(sd_gen,'r','FaceAlpha',0.8);
% bar(sd_sim,'g','FaceAlpha',0.8); 
% hold off
% legend('HRIR generica', 'HRIR proposta', 'Location', 'best')
% title('Distorï¿½ï¿½o Espectral Logarï¿½tmica')
% xlabel('Azimute [ï¿½]');
% ylabel('LSD [dB]')
% %x tick labels
% ticknum = dtf_med.SourcePosition(idxmed,1);
% ticknum = ticknum(1:4:end);
% xticks(1:4:length(idxmed));
% xticklabels(cellstr(num2str(ticknum)));
% 
% %%% RMSE %%%
% figure()
% bar(rmse_sim,'g','FaceAlpha',0.8); 
% hold on 
% bar( rmse_gen,'r','FaceAlpha',0.8);
% hold off 
% legend('HRIR proposta', 'HRIR genï¿½rica', 'Location', 'best')
% title('Erro RMS da resposta impulsiva')
% xlabel('Azimute [ï¿½]');
% ylabel('RMSE [dB]')
% xticks(1:4:length(idxmed));
% xticklabels(cellstr(num2str(ticknum)));
% 
% clc
% % Mï¿½dia e desvio padrï¿½o do erro no plano horizontal
% media_RMSE_generic = mean(rmse_gen)
% dp_rmse_generic = sqrt(var(rmse_gen))
% % 
% media_RMSE_sim = mean(rmse_sim)
% dp_rmse_prop = sqrt(var(rmse_sim))
% 
% media_LSD_genaric = mean(sd_gen)
% media_LSD_sim = mean(sd_sim)
