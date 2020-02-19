% EXPLORAR REDUNDANCIA DE DAS HRTFs
% Considerando que para angulos simetricos em relação a 0°
% (ex.: -30° e 30°)
% a resposta observada para cada orelha será virtualmente a mesma.
% (ex.: a resposta vista pela orelha esquerda a um angulo de -30°
%        será a mesma vista pela orelha direita a um angulo de 30°).
%
% A função explora essa redundancia para aumentar o número de casos para
% posterior treinamento em redes neurais.

%%% INPUT %%%
% SOFA object  (PS.: make sure the coordinate system is spherical
%                                                azi: 0° -> 360°
%                                                ele: -90° -> 90°). 

% Obj = dtf_ITA(1).dados;
function [Obj1, Obj2] = hrtf_simetry_split(Obj)
%% Spherical 2 Navigational
[no_directions, ~, ~] = size(Obj.Data.IR);
for l = 1:no_directions   
    azi = Obj.SourcePosition(l, 1); 
    ele = Obj.SourcePosition(l, 2); 

    [azi,ele] = sph2nav(azi,ele);

    meta.L(l, 1) = azi; 
    meta.L(l, 2) = ele; 
end

%% split IR in 2 Obj
objdata.L = Obj.Data.IR(:,1,:); 
objdata.R = Obj.Data.IR(:,2,:);
meta.R = meta.L;
%% Find correspondent azimuth in the same elevation and exchange positions
% there will be only one symmetric position
meta.R(:,1) = (abs(meta.R(:,1)));
idxel = nan(length(meta.R), 2);
for k = 1:length(meta.R)
        % make sure is the same elevation  
          idxel(k, :) = (find((meta.R(:,2) == meta.R(k,2)) & ...
                               meta.R(:,1) == meta.R(k,1)));
end
idxel = unique(idxel, 'rows');

% swap positions only in the IR
for k = 1:length(idxel)
    pos = idxel(k,:);
    objdata.R([pos(1), pos(2)], :) = objdata.R([pos(2), pos(1)], :);
end

%% PREPARE OUTPUT
Obj1 = Obj; Obj2 = Obj;
Obj1.SourcePosition = Obj.SourcePosition;
Obj2.SourcePosition = Obj.SourcePosition;
Obj1.Data.IR = objdata.L; 
Obj2.Data.IR = objdata.R; 

%% Plot symmetri between channels

% type = 'MagHorizontal';
% subplot(211)
% SOFAplotHRTF(Obj1, type);
% subplot(212)
% SOFAplotHRTF(Obj2, type);

end