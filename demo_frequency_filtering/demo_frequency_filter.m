function [output_ims] = demo_frequency_filter(input,lowCutOff,highCutOff,sigma)
% [output_ims] = demo_frequency_filter(input,lowCutOff,highCutOff,sigma)
% input is the only required variable
%
% code modified from SFFiltering1script.m by Izumi Ohzawa
% https://visiome.neuroinf.jp/database/item/6106
%
% this code will demonstrate how low pass/high pass/band pass filtering
% works. this code will output figures of your input image filtered in
% different ways. i recommend going through the code to understand the math
% behind filtering in fourier space (something i'm still trying to
% understand...). 
%
% variables:
% 
%   input - an image file loaded in using imread.
%   
%   lowCutOff - lowest frequency to cut off in bandpass demo. in
%   cycles/pixel.
%       recommend 0.02 (default) (lowest possible value - 0)
%   
%   highCutOff - highest frequency to cut off in bandpass demo. in
%   cycles/pixel.
%       recommend 0.2 (default) (highest possible value - 0.5) 
%  
%   sigma - width of gaussian function
%       recommend 0.05 (default)
%
%   output_ims - all images produced by this script will be put in this
%   cell variable (may have diff centre SF and sigma based on input)
%       output_ims{1} - lowpass filter
%       output_ims{2} - highpass filter
%       output_ims{3} - bandpass, centre SF 0.02, gaussian sigma 0.05
%       output_ims{4} - bandpass, centre SF 0.0650, gaussian sigma 0.05
%       output_ims{5} - bandpass, centre SF 0.1100, gaussian sigma 0.05
%       output_ims{6} - bandpass, centre SF 0.1550, gaussian sigma 0.05
%       output_ims{7} - bandpass, centre SF 0.2000, gaussian sigma 0.05
%
% the units of SF here correspond to cycles/pixel. The largest value
% possible is 0.5 (highest frequency), because it takes 2 pixels to have at
% least 1 cycle (black to white).
%
% the range used for the bandpass filter is 0.02 to 0.20. higher values
% (eg 0.5) were not used because weird edge artefacts occur when applying a
% filter at the edge of the spectrum. to get around this, you would have to
% considering padding your image.
%
% dependent scripts:
%       -GetLuminanceImage.m (by Izumi Ohzawa)
%           https://visiome.neuroinf.jp/database/item/6106
%       -Myff2.m (by Izumi Ohzawa)
%           https://visiome.neuroinf.jp/database/item/6106
%
% log:
% 20200311 - initialised. zoeyisherwood.
% contact: zoey.isherwood@gmail.com

%% start processing...-----------------------------------------------------
% check inputs-------------------------------------------------------------

if ~exist('lowCutOff','var')
    lowCutOff = 0.02;
end

if ~exist('highCutOff','var')
    highCutOff = 0.2;
end

if ~exist('sigma','var')
    sigma = 0.05;
end

% check if rgb and if dims are a factor of two-----------------------------

[outputim,maskmap] = makeimagefactortwo(input);

[row_start,cols_start] = find(maskmap, 1, 'first');
[row_end,cols_end] = find(maskmap, 1, 'last');


% set parameters-----------------------------------------------------------

% calculate the number of points for FFT (power of 2)
FFT_pts = 2 .^ ceil(log2(size(outputim)));

[A,~,~,mfx,mfy] = Myff2(outputim, FFT_pts(1), FFT_pts(2));

SF = sqrt(mfx .^ 2 + mfy .^ 2);

% low pass and high pass filter first...-----------------------------------

% parameters for SF filtering
ctrSF = 0; %0.06; % center SF

% SF-bandpass and orientation-unselective filter
lowpass_filt = exp(-(SF - ctrSF) .^ 2 / (2 * sigma ^ 2));
highpass_filt = 1-lowpass_filt;

A_lowpass_filtered = lowpass_filt .* A; % SF filtering
L_lowpass_filtered = real(ifft2(ifftshift(A_lowpass_filtered))); % IFFT
L_lowpass_filtered = L_lowpass_filtered(1: size(L_lowpass_filtered, 1), 1: size(L_lowpass_filtered, 2));
L_lowpass_filtered = uint8(rescale(L_lowpass_filtered).*255);

% crop images back to original size if they were padded--------------------

L_lowpass_filtered = L_lowpass_filtered(row_start:row_end,cols_start:cols_end);
lowpass_filt = lowpass_filt(row_start:row_end,cols_start:cols_end);

output_ims{1} = L_lowpass_filtered;

A_highpass_filtered = highpass_filt .* A; % SF filtering
L_highpass_filtered = real(ifft2(ifftshift(A_highpass_filtered))); % IFFT
L_highpass_filtered = L_highpass_filtered(1: size(L_highpass_filtered, 1), 1: size(L_highpass_filtered, 2));
L_highpass_filtered = uint8(rescale(L_highpass_filtered).*255);

% crop images back to original size if they were padded--------------------

L_highpass_filtered = L_highpass_filtered(row_start:row_end,cols_start:cols_end);
highpass_filt = highpass_filt(row_start:row_end,cols_start:cols_end);
output_ims{2} = L_highpass_filtered;

% plot high/low pass images with their masks-------------------------------

figure;

% original image
subplot(2,3,1);imshow(input); title('input image')

% low pass filtered image
subplot(2,3,2);imshow(L_lowpass_filtered); title('low pass filtered')

% high pass filtered image
subplot(2,3,3);imshow(L_highpass_filtered);title('high pass filtered')

% amplitude spectrum of input

min_amp_to_show = 10 ^ -10; % small positive value to replace 0 for log SF spectrum
ampspectra = abs(A);

show_log_amp = true;

if show_log_amp
    ampspectra(find(ampspectra < min_amp_to_show)) = min_amp_to_show; % avoid taking log 0
    ampspectra = log10(ampspectra);
end

ampspectra = uint8(rescale(ampspectra).*255);

ampspectra = ampspectra(row_start:row_end,cols_start:cols_end);

subplot(2,3,4);imshow(ampspectra); title('amplitude spectrum');

% low pass filter
subplot(2,3,5);imshow(lowpass_filt); title('low pass gaussian filter')

% high pass filter
subplot(2,3,6);imshow(highpass_filt); title('high pass gaussian filter')

suptitle('Low and High pass filtering')

% now bandpass filter------------------------------------------------------

% parameters for SF filtering

nsteps = 5;
im_counter = 2;
counter = 0;

figure;

for ctrSF = linspace(lowCutOff,highCutOff,nsteps) %0.06; % center SF
    
    im_counter = im_counter + 1;
    counter = counter + 1;
    
    % SF-bandpass and orientation-unselective filter
    filt = exp(-(SF - ctrSF) .^ 2 / (2 * sigma ^ 2));
    
    A_filtered = filt .* A; % SF filtering
    L_filtered = real(ifft2(ifftshift(A_filtered))); % IFFT
    L_filtered = L_filtered(1: size(L_filtered, 1), 1: size(L_filtered, 2));
    L_filtered = uint8(rescale(L_filtered).*255);
    
    L_filtered = L_filtered(row_start:row_end,cols_start:cols_end);
    filt = filt(row_start:row_end,cols_start:cols_end);
    
    output_ims{im_counter} = L_filtered;
    
    subplot(2,nsteps,counter);imagesc(L_filtered);colormap(gray);
    title(['SF center: ' num2str(ctrSF) ' sigma: ' num2str(sigma)]);
    subplot(2,nsteps,counter+nsteps);imagesc(filt);colormap(gray);
    
end

suptitle('Bandpass filtering')


end