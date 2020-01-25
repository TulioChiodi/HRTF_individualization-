function Obj_out = fit_2_CIPIC_grid(Obj_in, varargin)
% Converte objeto sofa qualquer aos parametros de medição utilizados nas
% medições do dataset CIPIC
% 
%   Input Parameters:
%    Obj_in:     Objeto de HRTFs SOFA com coordenadas esféricas
%                azi:    0° -> 360°     
%                elev: -90° -> 90°     

%   Output Parameters:
%     Obj_out:   Objeto de HRTFs SOFA com as característica 
%                de medição do dataset CIPIC

%   'fit_2_CIPIC_grid' aceita os seguintes parametros:
%     
%    'interp':   interpola as posicoes de interesse com a funcao interna 
%                'interpolateHRTF', necessaria Audio Toolbox
%    'move':     seleciona as posicoes mais proximas do grid objetivo e
%                força a assumirem suas coordenadas

%% Parse Arguments
p = inputParser;
defaultMethod = 'move';
validMethods = {'move','interp'};
checkMethod = @(x) any(validatestring(x,validMethods));
    
addRequired(p,'Obj_in',@isstruct);
addOptional(p,'method',defaultMethod,checkMethod)
parse(p,Obj_in,varargin{:})

%% Sample rate match
fs_in = Obj_in.Data.SamplingRate;
Fs = 44100; % sample rate no dataset CIPIC
if fs_in ~= Fs
    Obj_in = sofaResample(Obj_in, Fs);
end

%% Definir Grid Objetivo %%
% Converter coordenadas horizontais (CIPIC style) para esféricas (SOFA)
lat1 = [-80 -65 -55 -45:5:45 55 65 80];    % lateral angles (azimutes)
pol1 = -45 + 5.625*(0:49);                % polar angles (elevações)
pol  = repmat(pol1',length(lat1),1);
ida  = round(0.5:1/length(pol1):length(lat1)+0.5-1/length(pol1));
lat  = lat1(ida);
ii=1;
for aa = 1:length(lat1)
	for ee = 1:length(pol1)
		[azi,ele] = hor2sph(lat(ii),pol(ii));
        cipic.SourcePosition(ii,:) = ([azi, ele, 1]);
		ii = ii+1;
	end
end

% Caso queira fazer fit para outro grid qualquer, defina 'pos' com as
% coordenadas (esfericas) do grid objetivo 
meta.pos = Obj_in.SourcePosition;
pos      = cipic.SourcePosition;


%% Find new grid positions (MOVE)
switch p.Results.method
    case validMethods{1}
        idx = zeros(length(pos), 1);
        meta.fittedPOS = zeros(size(pos));
        for zz = 1:length(pos)  
            % Find nearest position of CIPIC data in the INPUT grid (Obj.SourcePosition)
            tsqr = sqrt((meta.pos(:,1)-pos(zz,1)).^2 + (meta.pos(:,2)-pos(zz,2)).^2);
            [~,idx(zz)] = min(tsqr); 
            % 'é possivel que uma mesma posicao no grid original seja escolhida 
            %para mais de uma posicao no grid referencia'

            % Posicoes selecionadas no grid original (util para visualização)
        %     meta.fittedPOS(zz,:) = Obj_in.SourcePosition(idx(zz),:);
        end
        % Assumindo posição no grid objetivo
        % (caso for visualizar as posições selecionadas, comentar a linha abaixo)
        meta.fittedPOS = cipic.SourcePosition;
        meta.fittedIR  = Obj_in.Data.IR(idx, :, :);

        % sort data (facilita a analise comparativa)
        [~, sort_idx]  = sortrows(round(meta.fittedPOS), [2, 1], {'descend', 'ascend'});
        meta.fittedPOS = meta.fittedPOS(sort_idx,:);
        meta.fittedIR  = meta.fittedIR(sort_idx,:,:);

%% Interpolar posicoes
    case validMethods{2}
        meta.fittedIR = interpolateHRTF(Obj_in.Data.IR, meta.pos(:, 1:2), pos(:, 1:2),...
                                       'Algorithm', 'vbap');   
        meta.fittedPOS = pos;
end

%% OUTPUT data (assembly and metadata) 
Obj_out = Obj_in;
Obj_out.Data.IR = meta.fittedIR;
Obj_out.SourcePosition = meta.fittedPOS;
if isfield(Obj_out, 'MeasurementSourceAudioChannel')
    Obj_out = rmfield(Obj_out, 'MeasurementSourceAudioChannel');
end
if isfield(Obj_out, 'MeasurementAudioLatency')
Obj_out = rmfield(Obj_out, 'MeasurementAudioLatency');
end

warning('off','SOFA:upgrade');
Obj_out = SOFAupgradeConventions(Obj_out);
Obj_out = SOFAupdateDimensions(Obj_out);


%% Plots
% %%% plot input %%%
% Obj_in.SourceView = Obj_in.SourcePosition;
% Obj_in.SourceView_Type = 'spherical';
% Obj_in.API.Dimensions.SourceView  = 'MC';
% SOFAplotGeometry(Obj_in)
% view([35 20])
% legend off
% 
% % %%% plot output %%%
% Obj_out.SourceView = Obj_out.SourcePosition;
% Obj_out.SourceView_Type = 'spherical';
% Obj_out.API.Dimensions.SourceView  = 'MC';
% SOFAplotGeometry(Obj_out)
% view([35 20])
% legend off
end