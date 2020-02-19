function Obj_out = sofaFit2Grid(Obj_in, out_pos, varargin)
% Converte posições de HRIRs SOFA às posições especificadas em 'out_pos'
% Davi R. Carvalho @UFSM - Engenharia Acustica - fevereiro/2020

%   Input Parameters:
%    Obj_in:     Objeto de HRTFs SOFA com coordenadas esféricas 
%                azi:    0° -> 360°     
%                elev: -90° -> 90° 
%                (com o último update, talvez funcione com outros sistemas de coordenadas)
%    out_pos:    Nx3 matrix de posições desejadas, em que N corresponde ao 
%                número total de posições, e as colunas correspondem a azimute,
%                elevação e raio respectivamente 

%   Output Parameters:
%     Obj_out:   Objeto de HRTFs SOFA com as característica 
%                de medição do dataset CIPIC

%   'fit_2_CIPIC_grid.m' aceita os seguintes parametros opcionais:
%     
%    'move':     Seleciona as posicoes mais proximas do grid objetivo e
%                força a assumirem suas coordenadas (metodo Default).
%    'interp':   Interpola as posicoes de interesse com a funcao interna 
%                'interpolateHRTF', necessaria Audio Toolbox.
%    'Fs':       Taxa de amostragem no objeto e saída (Default: 44100).
%
% Matlab R2019a
%% Parse Arguments
defaultMethod = 'move';
validMethods = {'move','interp'};
checkMethod = @(x) any(validatestring(x,validMethods));

defaultInterpMethod = 'bilinear';
validInterpMethods = {'bilinear','vbap'};
checkInterpMethod = @(x) any(validatestring(x,validInterpMethods));

p = inputParser;
addRequired(p,'Obj_in',@isstruct);
addOptional(p,'method',defaultMethod,checkMethod)
addOptional(p,'interpMethod',defaultInterpMethod,checkInterpMethod)

paramName = 'Fs';
defaultVal = 44100;
addParameter(p,paramName,defaultVal)
parse(p,Obj_in,varargin{:})

%% Sample rate match
fs_in = Obj_in.Data.SamplingRate;
if fs_in ~= p.Results.Fs
    Obj_in = sofaResample(Obj_in, p.Results.Fs);
end

%% Find new grid positions (MOVE)
switch p.Results.method
    case validMethods{1}
        meta.pos = Obj_in.SourcePosition;
        idx = zeros(length(out_pos), 1);
        meta.fittedPOS = zeros(size(out_pos));
        for zz = 1:length(out_pos)  
            % Find nearest position of CIPIC data in the INPUT grid (Obj.SourcePosition)
            tsqr = sqrt((meta.pos(:,1)-out_pos(zz,1)).^2 + (meta.pos(:,2)-out_pos(zz,2)).^2);
            [~,idx(zz)] = min(tsqr); 
            % 'é possivel que uma mesma posicao no grid original seja escolhida 
            %para mais de uma posicao no grid referencia'

            % Posicoes selecionadas no grid original (util para visualização)
        %     meta.fittedPOS(zz,:) = Obj_in.SourcePosition(idx(zz),:);
        end
        % Assumindo posição no grid objetivo
        % (caso for visualizar as posições selecionadas, comentar a linha abaixo)
        meta.fittedPOS = out_pos;
        meta.fittedIR  = Obj_in.Data.IR(idx, :, :);

        % sort data (facilita a analise comparativa)
%         [~, sort_idx]  = sortrows(round(meta.fittedPOS), [2, 1], {'descend', 'ascend'});
%         meta.fittedPOS = meta.fittedPOS(sort_idx,:);
%         meta.fittedIR  = meta.fittedIR(sort_idx,:,:);

%% Interpolar posicoes (INTERP)
    case validMethods{2}
        meta.pos = Obj_in.SourcePosition;
        meta.fittedIR = interpolateHRTF(Obj_in.Data.IR, meta.pos(:, 1:2), out_pos(:, 1:2),...
                                       'Algorithm', p.Results.interpMethod);   
        meta.fittedPOS = out_pos;
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