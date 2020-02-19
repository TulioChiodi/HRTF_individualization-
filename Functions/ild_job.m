function [hrir_l, hrir_r] = ild_job(hrir_l, hrir_r, xi)
% Estima ILD a partir de dados antropom�tricos
%%% INPUT %%% 
%% 
M = length(hrir_l); %n�mero de dire��es
%% Normliza��o das respostas impulsivas 
hrir_l = normalize(hrir_l, 'range');
hrir_r = normalize(hrir_r, 'range');

%% Calculo do ILD
D = 1:10; % ordem do modelo
a = []; %coeficientes
C = sum(a*xi);

for k = 1:n % frenquencia
    ild(k) = C(k)*sin(2*pi*k*m/M);
end
%% Aplicando diferen�a de n�vel
% if azi >= 0 
%     IR_l = hrir_l;
%     IR_r = hrir_r + ILD;
% else 
%     IR_l = hrir_l + ILD;
%     IR_r = hrir_r;
% end
end