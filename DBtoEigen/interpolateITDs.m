function itd = interpolateITDs(itds,ws)
% interpolateITDs linearly interpolates the ITDs specified in the vector 
% itds with the weights ws. 
%  
% SYNOPSIS: itd = interpolateITDs(itds,ws) 
% 
% INPUT itds: a vector of ITDs of size n 
%       ws: a vector of weights of length n used for the interpolation
%
% OUTPUT h: the resulting ITD  
% 
% REMARKS  
% 
% SEE ALSO 
% 
% AUTHOR    : Julian Villegas 
% $DATE     : 23-Mar-2017 14:00:24 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.1.0.441655 (R2016b) 
% FILENAME  : interpolateHRTFs.m 
%
n = size(itds,1);

for i=1:n
    itds(i,:) = itds(i,:) .* ws(i);
end
itd = sum(itds);
end 
% ===== EOF ====== [interpolateHRTFs.m] ======  
