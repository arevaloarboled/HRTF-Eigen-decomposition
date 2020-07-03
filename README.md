# HRTF-Eigen-decomposition
<p  align="center">
  <a  href="https://www.mathworks.com/products/matlab.html"  target="_blank">
    <img  src="https://img.shields.io/badge/Matlab-R2020a-blue.svg"  alt="Matlab R2020a">
  </a>
</p>

Audio spatialization in Matlab based on reconstruction of HRTFs from a Eigen decomposition. Synthetic HRTFs in the compressed database were set to have less than 1 dB spectral distortion between 0.1 and 16 kHz. The differences between the compressed measurements with those in the original database (Qu et al. databse [1]) do not seem to translate into degradation of perceptual location accuracy. 

## Usage
Input parameters:
* **Audio:** A monophonic audio vector sampled at 48 kHz.
* **Positions:**  A vector with the location path in spherical coordinates specifying the second in which the sound source will have the described location, the first position is taken as an initial position (distance, elevation, azimuth, seconds).
* **Output name (optional):** Filename (in `wav` format) to storage the binaural audio file resulting from the spatialization, by default, will be "spatialization.wav".
* **pinna (optional):** Specify the size of the pinna between Large 'l' or Small 's', by default is taken the Large pinna.

## Example

```matlab
%load an audio vector
[audio,fs]=audioread('monophonic_audio.wav');
%if the audio is not sampled at 48000, re sample to 48000
if fs!=48000
	audio=resample(audio,48000,fs);
end
%create the location path of the sound source
%this path goes through the horizon around the listener with distance 20, elevation 0
pos=[repelem(20,360);repelem(0,360);0:1:359;([1:size(audio,1)/360:size(audio,1)]/48000)]';
HRTF_Eigen_spatialization(audio,pos);
```

## Reference
  [1] T.  Qu,  Z.  Xiao,  M.  Gong,  Y.  Huang,  X.  Li,  and  X.  Wu,  “Distance-Dependent Head-Related Transfer Functions Measured With High Spa-tial  Resolution  Using  a  Spark  Gap,”IEEE Trans. on Audio, Speech &Language Processing, vol. 17, no. 6, pp. 1124–1132, 2009.
