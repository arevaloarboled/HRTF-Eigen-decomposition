function [sd,e] = sptrl_distor(s1,s2,n,fs)
%SPTRL_DISTOR Summary of this function goes here
%   Detailed explanation goes here
binRes = fs/n;
minFreq = 100;
maxFreq = 16000;
minBin = ceil(minFreq/binRes);
maxBin = floor(maxFreq/binRes);
if ~isreal(s1(1))
    mag1 = 20*log10(abs(s1));
    mag2 = 20*log10(abs(s2));
else
    mag1 = 20*log10(abs(fft(s1)));
    mag2 =  20*log10(abs(fft(s2)));
end
mag1 = mag1(minBin:maxBin);
mag2 = mag2(minBin:maxBin);
e = mag2(:)-mag1(:);
sd = (mean(e.^2))^0.5; % SD
end

