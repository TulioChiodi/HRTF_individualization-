function [IR_minL, IR_minR] = phase_job(hL, hR, itd, azi)
%%% input %%%%
% hL, hR: magnitudes para orelha esquerda e direita (espelhada e linear)
% itd: diferença iteraural de tempo (opcional)
% azi: azimute referente ao determinado itd (opcional)
N = length(hL);
half_hL = hL(1:N);
half_hR = hR(1:N);

%% Window spectrum 
retwin = ones(N*0.6, 1); 
hanwin = hanning((2*N)*0.4); % 10% de decaimento
win = [retwin; 
       hanwin(length(hanwin)/2:end-1)];
hL = half_hL.*win;
hR = half_hR.*win;

                %%% PHASE RECONSTRUCTION %%%
%% Minimum Phase 
phi_minL = imag(hilbert(-log(abs(hL))));
phi_minR = imag(hilbert(-log(abs(hR)))); 
% HRTF complexa
HminL = hL.*exp(1j*deg2rad(phi_minL));
HminR = hR.*exp(1j*deg2rad(phi_minR));

%% Back to time domain
IR_minL = real(ifft(HminL, N*2, 'symmetric'));
IR_minR = real(ifft(HminR, N*2, 'symmetric'));


%% Excess phase (ITD)
if nargin <= 2
    itd = 0;
    azi = 1;
end

offset = 10; 
if azi >= 0 
    IR_minL = circshift(IR_minL, round(itd+offset));
    IR_minL(1:round(itd+offset)) = 0;
    IR_minR = circshift(IR_minR, offset);
    IR_minR(1:offset) = 0;
else
    IR_minL = circshift(IR_minL, offset);
    IR_minL(1:offset) = 0;
    IR_minR = circshift(IR_minR, round(itd+offset));
    IR_minR(1:round(itd+offset)) = 0;
end



%% Metodo alternativo minimum phase
% for n = 1:2
%     if n == 1
%         hmag = hL;
%     else
%         hmag = hR;       
%     end
%     ctfflog=mean(log(abs(hmag)+eps),2);
%     Nfft = 2*N;
%     ctfcep = ifft(ctfflog, Nfft);
%     ctfcep(Nfft/2+2:Nfft) = 0;                % flip acausal part to causal part or simply multiply
%     ctfcep(2:Nfft/2) = 2*ctfcep(2:Nfft/2);    % causal part by 2 (due to symmetry)
%     ctfflog = fft(ctfcep,Nfft/2);
%     ctfp = exp(ctfflog);
% 
%     % get complex spectrum
% %     ctff = hmag .*exp(1j*angle(ctfp));
%     ctff = hmag .*exp(1j*deg2rad(ctfp));
% 
%     % get IR
%     IR.(['data', num2str(n)]) = (ifft(ctff, Nfft));
% end
% IR_minL = real(IR.data1);
% IR_minR = real(IR.data2);
end