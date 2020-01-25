function [itd, naz, nel] = itd_metodo3(wid, dep)
%%% Metodo 3 : ITD 
% Modelo geom�trico que simplifica a cabe�a a uma esfera 

%wid = distancia entre as duas orelhas;
% len = distancia entre face e nuca;

% REFERENCIA: Estimation of a Spherical-Head Model from Anthropometry 

% raio da esfera determinado a partir de uma soma ponderada das dimens�es da cabe�a 
% os pesos foram determinados a partir de uma regress�o sobre dados de
% varios indiv�duos

% Azimutes
naz = [-80, -65, -55 -45:5:45, 55, 65, 80]; 
% Eleva��es
nel = -45 + 5.625*(0:49); 

% HEAD DIMENSIONS 
wid = (wid/2);
dep = (dep/2);

% pesos
w1 = 0.51;
w2 = 0.18;
b = 3.2; 

% novo raio da esfera []
radius = (w1*wid + w2*dep + b);
c0 = 343e2; %[m/s]
% caso as dimensoes sejam desconhecidas recomenda-se a media da popula��o
% radius = 8.5e-2; % [m]

dim1 = length(naz);
dim2 = length(nel);
itd = zeros(dim1, dim2);
for k = 1:dim1
    for l = 1:dim2
        [theta, ~] = hor2sph(naz(k), nel(k)); 
        theta = deg2rad(theta);
        r = theta + sin(theta);
        itd(k, l) = (radius/c0)*r;
    end
end

fs = 44100;
itd(1:ceil(dim1/2), :) = flipud(itd(ceil(dim1/2):end,:));
itd = round(itd*fs);

%% PLOT
% figure()
% % surf(nel, naz, itd)
% surf(itd/fs*1e6)
% %%
% xlabel('eleva��o')
% xticks(1:5:length(nel))
% xticklabels(round(nel(1:5:end)))
% 
% ylabel('azimute')
% yticks(1:5:length(naz))
% yticklabels(naz(1:5:end))
% zlabel('tempo [\musec]')
% title('ITD - Indiv�duo (simulado)')
% view([0 0])
% 
% 
% %% medido
% % surf(ITD_norm(:,:,36)/fs*1e6)
% % xlabel('azimute')
% xticks(1:4:50)
% xticklabels(round(naz(1:2:end)))
% 
% ylabel('eleva��o')
% yticks(1:5:length(nel))
% yticklabels(round(nel(1:5:end)))
% zlabel('tempo [\musec]')
% title('ITD - Indiv�duo (medido)')
% view([0 0])
end


