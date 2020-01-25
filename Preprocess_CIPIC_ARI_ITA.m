clear all; clc; tic
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; Novembro/2019
%% General INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~INPUT SOFA FILE FROM THE CIPIC ARI AND ITA DATASETS,
% ~PUT THEM IN THE CIPIC GRID FORMAT AND SAME COORDINATE SYSTEM AND SAMPLE RATE
% ~APPLY SOFA DTF FUNCTION
%% PATHs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
addpath([pwd, '\Functions']);
% addpath(genpath('B:\Documentos\MATLAB\SOFA-master\API_MO')); %% Path com o toolbox SOFA

% ITA
pathita = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco ITA\SOFA\*.sofa');

% ARI
pathari = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco ARI\HRTF_DTF\database\ari\hrtf b_nh*.sofa');

% CIPIC
pathcipic = dir('B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\SOFA\*.sofa');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD CIPIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1 : length( pathcipic )
        CIPIC_ok(jj,1).dados = SOFAload( [pathcipic(jj).folder '\' pathcipic(jj).name], 'nochecks');
        CIPIC_ok(jj,1).dados = fit_2_CIPIC_grid(CIPIC_ok(jj,1).dados);
end
remove = [1:3, 5:8, 10, 36, 42]; %sem antropometria e (1 e 36) pra testes
CIPIC_ok(remove) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD ITA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
no_subj_ita = length(pathita);
for k = 1:no_subj_ita 
    ITA(k).dados = SOFAload([pathita(k).folder, '\',pathita(k).name], 'nochecks');
end

%%% Transição coordenadas cartesiana para esferica
for k = 1:no_subj_ita
    for l = 1:length(ITA(k).dados.SourcePosition)
        x = ITA(k).dados.SourcePosition(l, 1);  
        y = ITA(k).dados.SourcePosition(l, 2); 
        z = ITA(k).dados.SourcePosition(l, 3);
        % new coordinates
        [az,elev,r] = cart2sph(x,y,z);
        azi=rad2deg(az); elev=rad2deg(elev);
        [azi,ele]  = nav2sph(azi,elev);
        % update coordinates
        ITA(k).dados.SourcePosition(l, 1) = azi;
        ITA(k).dados.SourcePosition(l, 2) = ele; 
        ITA(k).dados.SourcePosition(l, 3) = round(r);
        % more metadata
        ITA(k).dados.SourcePosition_Type = 'spherical';
        ITA(k).dados.SourcePosition_Units = 'degree, degree, meter';              
    end
    % Adaptar ao grid CIPIC
    ITA_ok(k).dados = fit_2_CIPIC_grid(ITA(k).dados, 'interp');
end
% Remover indivíduo com dados incompletos
ITA_ok(14:15) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD ARI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
anthro_ari = load('anthro_ari.mat');
no_subj_ari = length(pathari); % TODAS AS HRTFS DISPONIBILIZADAS
for k = 1:no_subj_ari 
    id_ari_all(k) = sscanf(pathari(k).name, 'hrtf b_nh%d');  
end

% Encontrar hrtfs com antropometria
no_subj_ari = length(anthro_ari.id);
for k = 1:no_subj_ari 
    idx_import_ari(k) = find((anthro_ari.id(k) - 3000) == id_ari_all);
end

% Carregar HRTFs
for k = 1:no_subj_ari 
    idx = idx_import_ari(k);
    ARI(k).dados = SOFAload([pathari(idx).folder, '\',pathari(idx).name], 'nochecks');
end

% Remover individuos com antropometria INCOMPLETA
load('remove_ari.mat');
ARI(col_loc) = [];

% Reorganizar coordenadas e Resample (CIPIC style)
for k = 1:length(ARI)
    ARI_ok(k).dados = fit_2_CIPIC_grid(ARI(k).dados, 'interp');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HRTF 2 DTF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[no_directions, no_channels, no_samples] = size(CIPIC_ok(1).dados.Data.IR);
tic

% ITA
for k = 1:length(ITA_ok)
    [dtf_ITA(k).dados, ~] = SOFAhrtf2dtf(ITA_ok(k).dados, 'log');
     PTF_ITA(:,:,:,k) = 20*log10(abs(fft(dtf_ITA(k).dados.Data.IR, no_samples, 3)));
end
 PTF_ITA = permute(PTF_ITA, [3,4,1,2]);
 
% ARI
for k = 1:length(ARI_ok)
    [dtf_ARI(k).dados, ~] = SOFAhrtf2dtf(ARI_ok(k).dados, 'log');

    PTF_ARI(:,:,:, k) = 20*log10(abs(fft(dtf_ARI(k).dados.Data.IR, no_samples, 3)));
end
  PTF_ARI = permute(PTF_ARI, [3,4,1,2]);

% CIPIC
for k = 1:length(CIPIC_ok)
    [dtf_CIPIC(k).dados, ~] = SOFAhrtf2dtf(CIPIC_ok(k).dados,'log');

    PTF_CIPIC(:,:,:, k) = 20*log10(abs(fft(dtf_CIPIC(k).dados.Data.IR, no_samples, 3)));
end
 PTF_CIPIC = permute( PTF_CIPIC, [3,4,1,2]);

PTF = (cat(2, [PTF_CIPIC, PTF_ARI, PTF_ITA]))+100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero phase filter
N = no_samples/2;
retwin = ones(N*0.7, 1); 
hanwin = hanning((2*N)*0.3); % 10% de decaimento
win1 = [retwin;
        hanwin(length(hanwin)/2:end-1)];
win2 = flipud(win1);
win = [win1; win2];
% plot(win);

% Apply filter
PTF = bsxfun(@times,PTF, win);

%Normalize
PTF = PTF - 80;
surf(PTF(:,:,1,1), 'Linestyle', 'none')

% SAVE
save('PTF_cipic_ari_ita.mat', 'PTF')
disp('Dados Salvos!')


%% Export HRIRs: treinamento deep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unir bancos
% for jj = 1 : length( CIPIC_ok )
%     hrir(jj).dados = CIPIC_ok(jj).dados.Data.IR;  
% end
% for k = 1 : length( ARI_ok )
%     hrir(jj + k).dados = (ARI_ok(k).dados.Data.IR(:,:, 1:200))*20;  
% end
% 
% 
% dim1 = 25; dim2 = 50; no_channels = 2;
% for n = 1:length(hrir)
%     for i = 1:no_channels
%         cont = 0;
%         for k = 1:dim1
%             for l = 1:dim2
%                 cont = cont+1;
%                 target(k,l, :, n, i) = hrir(n).dados(cont, i, :);
%             end
%         end
%     end
% end

% save('target_deepnet.mat', 'target');
%%
% % PLOT ilustração da diferença entre hrtfs para diferentes pessoas numa
% mesma posição
%
%
% fs = 44100;
% Nt = 200;  
% % f_sw = linspace(0,(Nt-1)*fs/Nt, Nt);
% f_sw = linspace(0, fs-fs/Nt, Nt);
% h1 = 20*log10(abs(fft(squeeze(CIPIC_ok(1,1).dados.Data.IR(609, 2, :))))); % az = 0, ele = 0
% h2 = 20*log10(abs(fft(squeeze(CIPIC_ok(5,1).dados.Data.IR(609, 2, :))))); % az = 0, ele = 0
% h3 = 20*log10(abs(fft(squeeze(CIPIC_ok(11,1).dados.Data.IR(609, 2, :))))); % az = 0, ele = 0
% 
% %%%% PLOT %%%%
% figure()
% plot(f_sw(1:Nt/2), h1(1:Nt/2), 'lineWidth', 1.7, 'color',[0,0,0]); hold on
% plot(f_sw(1:Nt/2), h2(1:Nt/2), 'lineWidth', 1.7, 'color',[0,0.4,1]); 
% plot(f_sw(1:Nt/2), h3(1:Nt/2), 'lineWidth', 1.7, 'color',[1,0,0]); hold off
% 
% legend('ID 003','ID 011', 'ID 020', 'Location', 'best',...
%                                         'lineWidth', 1.5)
%     
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [dB]');
% axis([230 2e4 -30 20])
% title('HRTFs lado direito para azimute 0° e elevação 0°');
% set(gca,'FontSize',13)
% 
% grid on 


%% Split Symmetry
% % for k = 1:length(dtf_ITA)
% %         [Obj1, Obj2] = hrtf_simetry_split(dtf_ITA(k).dados);
% %         dtf_ITA_ok(k,1).dados = Obj1;
% %         dtf_ITA_ok(k,2).dados = Obj2;
% % end
% % 
% % for k = 1:length(dtf_ARI)
% %         [Obj1, Obj2] = hrtf_simetry_split(dtf_ARI(k).dados);
% %         dtf_ARI_ok(k,1).dados = Obj1;
% %         dtf_ARI_ok(k,2).dados = Obj2;
% % end
% % 
% % for k = 1:length(dtf_CIPIC)
% %         [Obj1, Obj2] = hrtf_simetry_split(dtf_CIPIC(k).dados);
% %         dtf_CIPIC_ok(k,1).dados = Obj1;
% %         dtf_CIPIC_ok(k,2).dados = Obj2;
% % end 
