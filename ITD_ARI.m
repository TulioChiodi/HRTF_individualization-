clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
% Calculo do ITD para - ARI dataset

%% PARTE 1: Carregar o banco de HRIRs  

path = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas\ARI_SOFA_ITD_EVAL\';
addpath(genpath(path));
addpath(genpath('B:\Documentos\MATLAB\SOFA-master\API_MO')); %% Path com o toolbox SOFA


anthro = load('anthro_ari.mat');
[id, idx] = sort(anthro.id);
%% Carregando dataset 
for k = 1:60 
    data = [path, 'individuo_', num2str(k) '.sofa'];
    dataset_itd(k).dados = SOFAload(data, 'nochecks');
end

% volta a ordem que tem nos dados anthropométricos,
% para quando aplicar na RNA ter os dados compatíveis
dataset_itd = dataset_itd(idx);
[~, no_channels, no_samples] = size(dataset_itd(1).dados.Data.IR);
no_subjects = length(dataset_itd);


%% Extração da HRIR
% 
% for k = 1:no_subjects
%     hrir(k).IR = dataset_itd(k).dados.Data.IR;
% end
% fs = dataset_itd(k).dados.Data.SamplingRate;


%% Conversão de SOFA >> CIPIC format
% conversão a fim de ajustar direções iguais ao do CIPIC
for n = 1:no_subjects
    % SOFA >> Formato ARI
    [hM, meta, stimPar] = SOFAconvertSOFA2ARI(dataset_itd(n).dados);

    % Formato ARI com direções CIPIC
    showfigures = 0;
    [hrir_l, hrir_r, name] = AA_ARI2CIPIC(hM, meta, stimPar, showfigures);    
    data_ari_cipic(n).hrir_l = hrir_l;
    data_ari_cipic(n).hrir_r = hrir_r;
    data_ari_cipic(n).name = name;
    
end
[dim1, dim2, no_samples] = size(data_ari_cipic(1).hrir_r);


%% ITD metodo ITA_start_IR
A = itaAudio; B = itaAudio;
% ITD_ita = zeros(no_directions, no_subjects);

for m = 1:no_subjects
    for k = 1:dim1
        for l = 1:dim2
            A.time =  squeeze((data_ari_cipic(m).hrir_l(k,l,:)));
            B.time =  squeeze((data_ari_cipic(m).hrir_r(k,l,:)));

            IL = ita_start_IR(A);
            IR = ita_start_IR(B);
        
%             temp = 1e6*abs(IL - IR)/fs; %valores em milissegundos 
            temp = abs(IL - IR);
            ITD_ita(k,l,m) = temp;
        end
    end
end


%% Corrigindo erros no ITD e smooth 

for k = 1:no_subjects
    for l = 1:dim2
        p = ITD_ita(13, l, k);
        temp = abs(ITD_ita(:, l, k) - p);
        ITDita_norm(:,l,k) = (temp);
%         ITD_ita = smooth(ITD_ita(:,dim2,k));
    end
end

%%
% ITD_ari = smoothdata(ITD_ita, 2);
ITD_ari = (ITDita_norm);

%% Remover indivíduos sem dados antropométricos 
load('remove_ari.mat')
ITD_ari(:,:, col_loc) =  [];


%% RESHAPE
[~,~,no_subjects] = size(ITD_ari);

for jj = 1:no_subjects
    for nt = 1:no_samples
        % esquerda
        H1 = ITD_ari(:,:,nt); % H1 é uma camada inteira da hrtf
        ITD_ari_f(nt, :, jj) = reshape(H1,[],1); %transformar a matriz em um vetor coluna
        % size(HRTF_d ) = [200 x 1250 x 36 x 2]
    end
end


%% SALVAR DADOS
save('DADOS_TREINAMENTO\ITD_ARI.mat','ITD_ari_f');
clc; disp('Dados salvos!')




