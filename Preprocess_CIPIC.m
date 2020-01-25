clear all; clc; 
% DAVI ROCHA CARVALHO; ENG. ACUSTICA - UFSM; MAIO/2019
% Analise antropométrica inicial, definição de parâmetros input para ANN
% por meio de Principal Component Analysis (PCA) - CIPIC dataset
cd 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\Rotinas'
tic
%% PARTE 1: Carregar o banco de HRIRs 
% Defina como path o local onde estão as pastas com os arquivos .mat do
% banco de dados.
onde = 'B:\Documentos\#3 - TCC\EAC-TCC-Davi\HRTF-datasets\Banco CIPIC\CIPIC_hrtf_database\standard_hrir_database';
path = dir([onde '\subject_*\hrir_final.mat']);
% Importa todo o banco de indivíduos para 'data'
for jj = 1 : length( path )
        data(jj) = importdata( fullfile( path(jj).folder, path(jj).name) );
end

% Remove indivíduos sem dados antropométricos e 1 indivíduo pra teste MANUAL
% individuos teste: 1, 36
remove = [1:3, 5:8, 10, 36, 42];
fields = fieldnames(data); %salva nome dos fields do struct para reconstrução 
data  = struct2cell(data); %transforma em celula para remover indivíduos
data(:,:, remove) = [];    %remove indivíduos 
data = cell2struct(data, fields, 1); %volta para struct


%% PARTE 2: Extração de dados [HRIR]
% Veja o grid de dispsição das caixas
no_subjects = length(data);
for jj = 1:no_subjects
    hrir(jj).L = data(jj).hrir_l; 
    hrir(jj).R = data(jj).hrir_r;
    
    temp = data(jj).ITD; 
    %transforma a matriz em um vetor coluna, em que todos os azimutes 
    % para 1 elevação são postos em sequencia.
    ITD_cipic(:, jj) = reshape(temp,[], 1); 
end


%% HRTF EM LOG SCALE [HRTF_log]
% Aplica fft e transforma para escala log
no_samples = length(hrir(1).L);
for jj = 1: no_subjects
    for m = 1:25
        for n = 1:50
            hrtf(jj).logL(m,n,:) = 20*log10(abs(fft(hrir(jj).L(m,n,:), no_samples))./ no_samples); 
            hrtf(jj).logR(m,n,:) = 20*log10(abs(fft(hrir(jj).R(m,n,:), no_samples))./ no_samples);
        end
    end
end
[dim1, dim2, no_samples] = size(hrtf(1).logL);

%% [HRTF_dir_log]

% Fazer a média de TODAS AS HRTFS DE UM INDIVIDUO e em seguida subtraimos a
% média de cada hrtf a fim de remover componentes que independem da
% direção
% A média é conhecida também como CTF (comum transfer function)
% A magnitude subtraida da média é conhecida como DTF (directional transfer functions)
no_directions = dim1*dim2; %1250 direções
no_channels = 2;

% somatório de todas as direções para determinado indivíduo e orelha
% calculating directional mean
med_vec = zeros(no_samples, no_subjects, no_channels);

for i = 1:no_subjects
    temp1 = zeros(no_samples,1); 
    for j = 1:dim1
       for k = 1:dim2
            temp1 = temp1 + squeeze(hrtf(i).logL(j, k, :)) + squeeze(hrtf(i).logR(j, k, :));
       end
    end
    temp1 = temp1/(2*no_directions);   
    med_vec(:, i) = temp1; % (no_samples x no_subjecs x no_channels)
    end

% "the mean is subtracted from each HRTFlog to obtain a new transfer
% function (HRTF_dir_log)
%  which represents primarily direction-dependent spectral effects"

for i = 1: no_subjects
    for j = 1:dim1
        for k = 1:dim2
            hrtf(i).dir_logL(j, k, :) = squeeze(hrtf(i).logL(j, k, :)) - med_vec(:, i);
            hrtf(i).dir_logR(j, k, :) = squeeze(hrtf(i).logR(j, k, :)) - med_vec(:, i);
        end
    end
end


%% Reestruturação dos dados para input na PCA
% Ref.: Dimensionality Reduction in HRTF by using Multiway Array Analysis

% Empilha cada camada 3d em um vetor coluna, 
% para cada indivíduo e direção teremos um vetor de samples
% naz = [azimute]
% nel = [elevação]
% nt  = [tempo]

%transforma a matriz em um vetor coluna, em que todos os azimutes 
% para 1 elevação são postos em sequencia.
% 
% for jj = 1:no_subjects
%     for nt = 1:no_samples
%         % esquerda
%         H1 = hrtf(jj).dir_logL(:,:,nt); % H1 é uma camada inteira da hrtf
%         HRTF_d(nt, :, jj, 1) = reshape(H1,[],1); %transformar a matriz em um vetor coluna
%         % size(HRTF_d ) = [200 x 1250 x 36 x 2]
%         
%         % direita
%         H2 = hrtf(jj).dir_logR(:,:,nt); 
%         HRTF_d(nt, :, jj, 2) = reshape(H2,[],1); %transformar a matriz em um vetor coluna     
%     end
% end

for jj = 1:no_subjects
    pos = 1;
    for k = 1:dim1
        for l = 1:dim2         
            HRTF_d(:, pos, jj, 1) = hrtf(jj).dir_logL(k, l, 1:no_samples/2);
            HRTF_d(:, pos, jj, 2) = hrtf(jj).dir_logR(k, l, 1:no_samples/2);
            pos = pos + 1;
        end
    end
end

%% Personal Transfer Function [PTF] 
%  realizar a média de HRTF_dir_log [HRTF_d] de todos os indivíduos para uma
% mesma direção (a fim de centralizar a matriz antes de aplicar a PCA)
%   Apesar de a funcao pca permitir centralizar automaticamente, a
%   centralizacao preia, permite a reobtenção dos dados originais
% med_vec2 = zeros(no_samples, no_directions, no_channels);
% for channel = 1:no_channels
%     for direction = 1:no_directions
%         temp3 = zeros(no_samples, 1);
%         for subject = 1:no_subjects
%              temp3 = temp3 + HRTF_d(:, direction, subject, channel);
%         end
%         temp3 = temp3/no_subjects;
%         med_vec2(:, direction, channel) = temp3;
%     end
% end
% 
% % Subtração da media nas HRTFs
% PTF_CIPIC = zeros(size(HRTF_d));
% for i = 1:no_channels
%     for j = 1:no_subjects
%         for k = 1:no_directions
%             PTF_CIPIC(:, k, j, i) = HRTF_d(:, k, j, i) - med_vec2(:, k, i);
%         end
%     end
% end
% A matrix precisa estar em determinada ordem para entrar na PCA
% "Rows of X correspond to observations and columns to variables"
PTF_CIPIC = permute(HRTF_d, [1, 3, 2, 4]); 

%% Smooth da PTF
% Josef recomenda aplicar smooth nos dados antes de aplicar a PCA, melhora
% a compressão nos dados por reduzir a variabilidade nos dados sem diminuir
% a precisão na localização
for k = 1:no_channels
    for l = 1:no_directions
        for m = 1:no_subjects
            PTF_CIPIC(:,m,l,k) = smooth(PTF_CIPIC(:,m,l,k), 'sgolay', 3);
        end
    end
end

 %% Principal Component Analysis (PCA) - METODO 1
no_PC = 10; %número de principais componentes de interesse
% [weight_vectors, PC_mtx2, eigen_value] = PCA_fun(PTF, no_PC);
% 12x 30x 1250x 2 - formato de input na neural net 

% Para calcular o quanto representa o numero de princcipais componentes
% representa o dataset inteiro, assumimos que com todos as PCs teremos 100%
% logo para normalizar é necessário somar todos os valores  de variancia e dividir o
% vetor original pela soma, assim ao somarmos os n primeiros valores,
% teremos a porcentagem do quanto aquela quantidade de PC representa o dataset original

% A PCA é aplicada para todos os sujeitos em uma mesma direção
for m = 1:no_channels
    for n = 1:no_directions        
        med_vec2(:,n,m) = mean(PTF_CIPIC(:, :, n, m),2);
        data_mtx = PTF_CIPIC(:, :, n, m) - med_vec2(:,n,m); % média subtraida de cada COLUNA  
         
        % Principal component analysis
        [coeff, score, latent] = pca(data_mtx, 'NumComponents', no_PC, ...
                                                  'Centered', false,...
                                                  'Algorithm','svd');
        %  "Component scores are the transformed variables
        % and loadings(coeff) describe the weights to multiply the normed original variable
        % to get the component score."
        % In this case pca does not center the data.
        % You can reconstruct the original data using score*coeff'       
        % Salva os dados de cada direção numa camada
        % Reconstrução e dada por (PC_mtx*weights)       
        weights(:,:, n, m) = coeff'; %Cada matrix coeff tem que ser orthonormal
        PC_mtx(:,:, n, m) = score;
        variancia(:, n, m) = latent;
    end
end
disp('PCA calculada!')

%% Dados ANTHROPOMETRICOS
% Carregar dados 
anthro = load('anthro.mat');
remove = [1:3, 5:8, 10, 36, 42];

% separando os 8 parametros necessários para a ann
% d1, d3 : d6
% d1 = anthro.D(:, [ 3, 4, 5, 6]);   % L
d1 = anthro.D(:, [1, 3, 4, 5, 6]);   % L
d2 = anthro.D(:, ([1, 3, 4, 5, 6] + 8)); % R
d = cat(3, d1, d2); %Concatena na terceira dimensão

% x1, x3, x12
x = anthro.X(:, [1, 3, 12]);

%removeremos os individuos que não foram medidos e indivíduos para teste
% indivíduos para teste: 1, 4, 9, 11:14.
d(remove , :, :) = []; 
x(remove , :) = []; 

% unir as duas matrizes
InputMatrix_CIPIC(:,:,1) = [d(:,:,1), x]';  % L
InputMatrix_CIPIC(:,:,2) = [d(:,:,2), x]';  % R

%% SALVAR DADOS
save('DADOS_TREINAMENTO\input_target_CIPIC.mat','InputMatrix_CIPIC', 'med_vec', 'med_vec2');
save('DADOS_TREINAMENTO\input_target.mat','InputMatrix_CIPIC', 'PC_mtx', 'weights', 'med_vec', 'med_vec2');

save('DADOS_TREINAMENTO\preprocess_CIPIC.mat', 'PTF_CIPIC', 'med_vec2');
disp('Dados salvos!')
toc

    
%% PLOT DO NUMERO DE PC NECESSÀRIOS 


% Para calcular o quanto representa o numero de princcipais componentes
% representa o dataset inteiro, assumimos que com todos as PCs teremos 100%
% logo para normalizar é necessário somar todos os valores  de variancia e dividir o
% vetor original pela soma, assim ao somarmos os n primeiros valores,
% teremos a porcentagem do quanto aquela quantidade de PC representa o dataset original

% numero de PC vs variancia
vari = variancia(:,1,1);
pn_pc = 1:length(vari(:,1,1));
var_norm = vari./sum(vari);
for i = 1:length(vari)
    var_var(i) = sum(var_norm(1:i));
end

figure(3)
plot(pn_pc, var_var*100, 'LineWidth', 2.0)
axis tight
ylabel('Representação do dataset original [%]')
xlabel('Numero de Componentes Principais [~]')
grid on 
% xtick = var_var*100;
% xticks 

