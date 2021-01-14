function hmin = minPhaseHRIR(h)
% MINPHASEHRIR converts an HRIR in its minimum-phase filter version 
%  
% SYNOPSIS: hmin = minPhaseHRIR(h)  
% 
% INPUT hrirs: a mxnx2  matrix of n hrirs with m taps 
% 
% OUTPUT hmin: the mxnx2 matrix in minimum-phase representation. 
% 
% REMARKS This could also be achieved using cepstral components. Maybe, in
% the future, we could investigate differences between the two methods.
% 
% SEE ALSO 
% 
% AUTHOR    : Julian Villegas 
% $DATE     : 24-Mar-2017 10:49:20 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.1.0.441655 (R2016b) 
% FILENAME  : minPhaseHRIR.m 
    hmin =  real(ifft(exp(conj(hilbert(log(abs(fft(h))))))));
end 
% ===== EOF ====== [minPhaseHRIR.m] ======  
