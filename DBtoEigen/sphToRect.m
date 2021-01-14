function p = sphToRect(aziEleDist)
% SPHTORECT translate spherical coordinates in degrees to rectangular 
% coordinates, approaching the result to the nearest integer. Effectively, 
% yielding a precision of 1cm in each dimension.
%  
% SYNOPSIS: p = sphToRect(aziEleDist)  
% 
% INPUT aziEleDist: an array of spherical coordinates (azimuth, elevation, 
% distance)  
% 
% OUTPUT p: the same array in rectangular coordinates: 
%       - positive x axis: front
%       - positive y axis: right
%       - positive z axis: above
% 
% REMARKS Azimuth in Qu's database is given clockwise  
% 
% SEE ALSO 
% 
% AUTHOR    : Julian Villegas 
% $DATE     : 13-Mar-2017 16:30:35 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.1.0.441655 (R2016b) 
% FILENAME  : sphToRect.m 
 [x,y,z] = sph2cart(deg2rad(aziEleDist(:,1)),... 
                               deg2rad(aziEleDist(:,2)),...
                               aziEleDist(:,3));
   %p = round([x,y,z]);    
   p = [x,y,z];    
end 
% ===== EOF ====== [sphToRect.m] ======  
