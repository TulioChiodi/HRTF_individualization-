function [xinpt, sig, mu] = nn_input_preprocess(InptMtx, sig, mu)
%%% Regularização de matriz entre 0 e 1 para input na ANN
% normalização dos dados de entrada a partir da média e da variância das
% features entre indivíduos 

%%% REFERENCIA %%% 
%Deep Neural Network Based HRTF Personalization Using
%Anthropometric Measurements
if nargin < 2  
    sig = std(InptMtx, 0, 2); %desvio padrão para cada parâmetro de entrada entre indivíduos
    mu  = mean(InptMtx, 2);   %média de cada parâmetro
end

xinpt = zeros(size(InptMtx)); % inicialização da matriz
[~, ~, no_channels] = size(InptMtx);
    for k = 1:no_channels
        xinpt(:,:,k) = (1 + exp(-(InptMtx(:,:,k) - mu(:,:,k))./sig(:,:,k))).^(-1);
    end
end