function [hrir_l, hrir_r] = phase_itd(IR_minL, IR_minR, itd )
%% Excess phase (ITD)
if nargin <= 2
    itd = 0;
    azi = 1;
end

%%
% Azimutes
naz = [-80, -65, -55 -45:5:45, 55, 65, 80]; 
% Elevações
nel = -45 + 5.625*(0:49); 

offset = 0;
% 
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
% isminphase(IR_minL)