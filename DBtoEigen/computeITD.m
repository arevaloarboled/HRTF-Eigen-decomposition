function [itd,tl,tr] = computeITD(h,sr)
% COMPUTEITD computes delays in seconds from the HRIRs in h.
%
% SYNOPSIS: [itd,tl,tr] = computeITD(h,sr)
%
% INPUT  h: a mxnx2 matrix of HRTFs where m is the number of FFT bins,
%           n is the number of HRTFs, the third dimension correspons to
%           the the left (1) and right channel (2) of a given HRTF.
%        sr: sampling frequency
%
% OUTPUT itd: array of inter-aural time delay. Positive values means that 
% the left channel is delayed 
%        tl: array of delays in s to the left ear
%        tr: array of delays in s to the right ear
%
%
% REMARKS Note the Threshold method is best according to my tests. This
% code is based on the one provided in:
% [1] J. Estrella, "On the extraction of interaural time differences from 
% binaural room impulse responses," Masters thesis, Technische Universitat
% Berlin, 2010.
%
% SEE ALSO testingITD.m
%
% AUTHOR    : Julian Villegas
% $DATE     : 23-Mar-2017 15:47:12 $
% $Revision : 1.00 $
% DEVELOPED : 9.1.0.441655 (R2016b)
% FILENAME  : computeITD.m
switch nargin
    case 0
        error('Insufficient arguments');
    case 1
        warning('no samplig rate defined. Assuming 2^16 Hz');
        sr = 2^16;
end

upfactor = 10;
tauUp = 1/(upfactor);
% -15 dB works fine for MIT dB
% -12 dB works better for Qu's database. See: testingITD.m
onset_threshold_dB = -12;
onset_threshold =10^(onset_threshold_dB/20);

upH = zeros(size(h,1)*upfactor,size(h,2),size(h,3));
for i=1:size(h,2)
    for j=1:size(h,3)
        upH(:,i,j) = upsample(h(:,i,j),upfactor);
    end
end

% create a matrix with max value, sample position, for each channel
maxVals = zeros(size(h,2),2,size(h,3));
for i=1:size(h,2)
    for j=1:size(h,3)
        [value, position] = max(abs(upH(:,i,j)));
        maxVals(i,:,j) = [value, position];
    end
end

results = zeros(size(h,2),3); % results as itd, left, and right delays

for i=1:size(h,2) % for each hrir
    for j=1:size(h,3) % for each channel
        % find the location where the tap is greater than the threshold 
        results(i,j+1) = find(abs(upH(1:maxVals(i,2,j),i,j)) > ...
            abs(maxVals(i,1,j)*onset_threshold),1);
    end
end

results(:,1) = results(:,2)-results(:,3);
results = results*tauUp/sr; % return samples in seconds

itd = results(:,1);
tl = results(:,2);
tr = results(:,3);
end
% ===== EOF ====== [computeITD.m] ======
