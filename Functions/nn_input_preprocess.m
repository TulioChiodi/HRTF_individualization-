function [xinpt, sig, mu] = nn_input_preprocess(InptMtx, sig, mu)
%%% Regulariza��o de matriz entre 0 e 1 para input na ANN
% normaliza��o dos dados de entrada a partir da m�dia e da vari�ncia das
% features entre indiv�duos 

%%% REFERENCIA %%% 
%Deep Neural Network Based HRTF Personalization Using
%Anthropometric Measurements
if nargin < 2  
    sig = std(InptMtx, 0, 2); %desvio padr�o para cada par�metro de entrada entre indiv�duos
    mu  = mean(InptMtx, 2);   %m�dia de cada par�metro
end

xinpt = zeros(size(InptMtx)); % inicializa��o da matriz
[~, ~, no_channels] = size(InptMtx);
    for k = 1:no_channels
        xinpt(:,:,k) = (1 + exp(-(InptMtx(:,:,k) - mu(:,:,k))./sig(:,:,k))).^(-1);
    end
end