function view_spectra(input)
% view_spectra(input)
%
% this code will produce figures of the following:
%   1) amplitude spectra,
%   2) phase spectra,
%   3) the amplitude spectra on a log log axis averaged across all
%   orientations
% 
% equation to combine magnitude and phase:
%   image_fourier = abs(amplitude_spectrum) .* exp(i * phase_spectrum);
%   image = ifft2(image_fourier);.
%
% log:
% 20200311: initialised. zji.
% contact: zoey.isherwood@gmail.com
%
%% start processing
% check if rgb and if dims are a factor of two-----------------------------

[outputim,maskmap] = makeimagefactortwo(input);

[row_start,cols_start] = find(maskmap, 1, 'first');
[row_end,cols_end] = find(maskmap, 1, 'last');

% get amp spectra and phase spectra----------------------------------------

spectrum = fftshift((fft2(double(outputim))));

amplitude_spectrum = abs(log10(spectrum));
% change values to fit between 0 and 255
amplitude_spectrum = uint8(rescale(amplitude_spectrum).*255);
%crop to original image size
amplitude_spectrum = amplitude_spectrum(row_start:row_end,cols_start:cols_end);

phase_spectrum = angle(spectrum);
% change values to fit between 0 and 255
phase_spectrum = uint8(rescale(phase_spectrum).*255);
%crop to original image size
phase_spectrum = phase_spectrum(row_start:row_end,cols_start:cols_end);

% plot everything----------------------------------------------------------

figure;

% input image----------------------------------------------------------

subplot(2,2,1);imshow(input);colormap(gray); title('input image');

% amplitude spectrum----------------------------------------------------
subplot(2,2,2); imshow(amplitude_spectrum);colormap(gray); title('amplitude spectrum')

% input image amp spectra log-log axis, averaged across orientations-------

subplot(2,2,3);
[f,s]=eo_polaraverage(abs(fft2(double(outputim))));
f2=log2(f);
s2=log2(s);
logf=f2(1):abs(f2(end)/length(f2)):f2(end);
logs=interp1(f2,s2,logf);
plot(logf,logs,'b','linewidth',3);
hold on;
xlabel('log f','FontSize',12,'FontName','Helvetica');
ylabel('log Amplitude','FontSize',12','FontName','Helvetica');
p_temp=polyfit(logf,logs,1);
plot(f2,p_temp(1)*f2+p_temp(2),'r','linewidth',3);
legend('image data','linear fitting');
set(gca,'FontSize',12,'FontName','Helvetica');
alpha=sprintf('alpha=%2.2f',p_temp(1));
text(logf(round(0.1*length(logf))),logs(round(0.8*length(logs))),alpha, ...
    'FontSize',12,'FontName','Helvetica');
axis equal

title('amplitude spectrum on log log axis')

% phase spectrum-----------------------------------------------
subplot(2,2,4); imshow(phase_spectrum);colormap(gray); title('phase spectrum')


end