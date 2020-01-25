clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
% Principal Component Analysis (PCA) 
tic
%% Carregando dados 
addpath('B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas\DADOS_TREINAMENTO\')
load('PTF_cipic_ari_ita.mat')
[no_samples, no_subjects, no_directions, no_channels]  = size(PTF);

%% CETER DATA 
% med_vec2 = squeeze(mean(PTFd, 2));
% 
% % Subtração da media nas HRTFs
% PTF = zeros(size(PTFd));
% for i = 1:no_channels
%     for j = 1:no_subjects
%         for k = 1:no_directions
%             PTF(:, j, k, i) = PTFd(:, j, k, i) - med_vec2(:, k, i);
%         end
%     end
% end

% %Smoothie  
% PTF = smoothdata(PTF, 1, 'sgolay');
surf(PTF(:,:,1,1), 'Linestyle', 'none')

%% Principal Component Analysis (PCA)
no_PC = 12; %número de principais componentes de interesse

% A PCA é aplicada para todos os sujeitos em uma mesma direção
for m = 1:no_channels
    for n = 1:no_directions
        med_vec2(:,n,m) = mean(PTF(:, :, n, m), 2);
        data_mtx = PTF(:, :, n, m) - med_vec2(:,n,m); % média subtraida de cada COLUNA  
         
        [coeff, score, variancia] = pca(data_mtx,'NumComponents', no_PC,  ...
                                                 'Centered', false,...
                                                 'Algorithm','svd');
        %  "Component scores are the transformed variables
        % and loadings(coeff) describe the weights to multiply the normed original variable
        % to get the component score."
        % In this case pca does not center the data.
        % You can reconstruct the original data using score*coeff'       
        % Salva os dados de cada direção numa camada
        weights(:,:, n, m) = coeff'; %Cada matrix coeff tem que ser orthonormal
        PC_mtx(:,:, n, m) = score;
    end
end
disp('PCA calculada!')


%% PLOT DO NUMERO DE PC NECESSÀRIOS 
% Para calcular o quanto representa o numero de princcipais componentes
% representa o dataset inteiro, assumimos que com todos as PCs teremos 100%
% logo para normalizar é necessário somar todos os valores  de variancia e dividir o
% vetor original pela soma, assim ao somarmos os n primeiros valores,
% teremos a porcentagem do quanto aquela quantidade de PC representa o dataset original

% numero de PC vs variancia 
pn_pc = 1:length(variancia);
var_norm = variancia./sum(variancia);
for i = 1:length(variancia)
    var_var(i) = sum(var_norm(1:i));
end

figure()
plot(pn_pc, var_var*100, 'LineWidth', 2.0)
axis tight
ylabel('Representação do dataset original [%]')
xlabel('Numero de Componentes Principais [~]')
grid on 
% 
% for k = 1:no_channels
%     for l = 1: no_directions
%         med_vec2(:,l,k) = (med_vec2(:,l, k) + med_vec2_ARI(:,l, k) + med_vec2_ITA(:,l, k))/3;
%     end
% end
%% SAVE data
save('DADOS_TREINAMENTO\target_cipic_ari_ita.mat', 'weights', 'PC_mtx', 'med_vec2')
disp('Dados salvos!')


%% VERIFY reconstruction 
% [no_samples, no_PC, no_directions, no_channels] = size(PC_mtx);

% for n = 1:no_channels
%     for i = 1:no_directions
%              recon(:,:, i, n) = PC_mtx(:, :, i, n)*squeeze((weights(:,:, i, n)));                                   
% %         hrtf_dir_log_sim(:, i, n) = hrtf_dir_log_sim(:, i, n) + med_vec2(:, i, n); 
%     end
% end
% 

% recon = PC_mtx*weights +med_vec2;


