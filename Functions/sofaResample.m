function [Obj] = sofaResample(Obj, Fs)

% sofaResample Resample (frequency) SOFA file
%
% Usage
%   [Obj] = sofaResample(Obj, Fs);
%
% Input
%   Obj : Sofa structure
%   Fs: target sampling frequency
%
% Output
%   Obj : resampled Sofa structure
%
% Authors
%   David Poirier-Quinot

% discard is already correct sampling frequency
if( Obj.Data.SamplingRate == Fs )
    warning('resampling discarded (already at correct sampling freq.)');
end

% resample
N = ceil( (Fs/Obj.Data.SamplingRate) * size(Obj.Data.IR, 3) ); % new length
IR = zeros(size(Obj.Data.IR, 1), size(Obj.Data.IR, 2), N);
for iPos = 1:size(Obj.Data.IR, 1)
    for iCh = 1:size(Obj.Data.IR, 2)
        IR(iPos, iCh, :) = resample(Obj.Data.IR(iPos, iCh, :), Fs, Obj.Data.SamplingRate);
    end
end
Obj.Data.IR = IR;

% update sampling rate
Obj.Data.SamplingRate = Fs;