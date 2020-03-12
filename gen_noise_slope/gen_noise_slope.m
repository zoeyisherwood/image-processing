function [output_noise] = gen_noise_slope(desired_slope,x_size,y_size,plotfigs)
% [output_noise] = gen_noise_slope(desired_slope,x_size,y_size,plotfigs)
% 
% variables:
%   
%   desired_slope - what you want the slope of the noise image to be
%       works best between 0 and -2.5. beyond this range, you may get
%       clipping artefacts.
% 
%   x_size - size of your desired noise image in the x direction
% 
%   y_size - size of your desired noise image in the y direction
% 
%   plotfigs - set to 1 if you want to see plots of the amplitude spectrum
%       before and after transformation. set to 0 if not.
% 
% this code will create a noise image with a specified amplitude slope
% value. any image dimension size/combination should work. however, note
% that if you use a non-square dimension, or if both dims are an odd number
% (eg 117x117), the noise image will be padded to be square/even before
% fourier transformation. after the slope of the noise has been manipulated
% and converted back to real space, the image will be unpadded. this is
% important to note since there may be slight discrepancies between the
% reported measured slope and the actual slope. this is due to the fact
% that the amplitude slope would be estimated on a slightly bigger image
% (due to padding). these discrepancies should be minimal but just
% keep it in mind. another reason the image is made square before fourier
% analysis is to ensure that the final image doesn't appear stretched
% out... the function that manipulates the slope needs equal dims.
% 
% % current limitations of script:
%   -only generates grayscale noise. plan to add compatibility
%   for colored noise...
%   -data only output in uint8 format, which can cause clipping artefacts
%   for steep slopes. plan to add higher bitrate options (2^10, 2^14) which
%   can be used on high bitrate monitors like Display++...
%   -x axis of amp spectra plots not very intuitive. change later...
%   -slope measurement is done in loglog space. want to eventually add an
%   option to measure the slope in linear space. this is particularly
%   important when rounding errors are emphasized when log transforming the
%   data and affect the fit. will do later...
%
% dependent scripts:
%   you need to have:
%       -makeimagesquare.m (contact zoeyisherwood for this)
%       -eo_polaraverage.m (part of one_over_f package)
%           https://visiome.neuroinf.jp/database/item/6110
% log:
% 20200311 - initialised by zoeyisherwood.
% 
% contact: zoey.isherwood@gmail.com

%% start processing:

%generate random noise-----------------------------------------------------
input = random('Normal',0,1,x_size,y_size); % generate image of Gaussian white noise

% make image square if it isn't already------------------------------------

% [squareim,maskmap] = makeimagesquare(input);
[squareim,maskmap] = makeimagefactortwo(input);

% get image size (it should be square)-------------------------------------

imsize = size(squareim,1);

% measure slope of input (needed to generate a slope of your choosing)-----

[f,s]=eo_polaraverage(abs(fft2(double(squareim))));
f2=log2(f);
s2=log2(s);
logf=f2(1):abs(f2(end)/length(f2)):f2(end);
logs=interp1(f2,s2,logf);
slope=polyfit(logf,logs,1);
slope=slope(1);

% change slope of input to slope of your choosing...-----------------------

[x,y]=meshgrid(linspace(-1,1,imsize));
r=sqrt(x.^2+y.^2);
q=r.^(-(-desired_slope+slope));
fftdata=fft2(double(squareim));
fftabs=abs(fftdata);
fftang=angle(fftdata);
m=ifftshift(q);
newfftabs=fftabs.*m;
m_data=real(ifft2(newfftabs.*exp(sqrt(-1)*fftang)));
m_data=m_data-min(m_data(:));
m_data=uint8(m_data/max(m_data(:))*255);

% measure slope of manipulated image---------------------------------------

[f,s]=eo_polaraverage(abs(fft2(double(m_data))));
f2=log2(f);
s2=log2(s);
logf=f2(1):abs(f2(end)/length(f2)):f2(end);
logs=interp1(f2,s2,logf);
manipulated_im_slope=polyfit(logf,logs,1);
manipulated_im_slope=manipulated_im_slope(1);

% crop images back to original size if they were padded--------------------

[row_start,cols_start] = find(maskmap, 1, 'first');
[row_end,cols_end] = find(maskmap, 1, 'last');

output_noise = m_data(row_start:row_end,cols_start:cols_end);

% report slopes------------------------------------------------------------

disp(['orig slope: ' num2str(slope)]);
disp(['manipulated slope: ' num2str(manipulated_im_slope)]);

% plot images + spectra if indicated---------------------------------------

if plotfigs == 1

    figure;
    
    % orig noise-----------------------------------------------------------
    
    if ndims(input) == 3
        input=input(:,:,1:3); %get rid of alpha channel if there.
    end
    
    subplot(2,2,1);imagesc(input);colormap(gray);
    title('input image');
    
    % manipulated noise----------------------------------------------------
    subplot(2,2,2); imagesc(output_noise);colormap(gray);
    title('manipulated image')
    
    % orig noise spectra---------------------------------------------------
    subplot(2,2,3);
    [f,s]=eo_polaraverage(abs(fft2(double(squareim))));
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
    alpha=sprintf('alpha=%2.2f',slope);
    text(logf(round(0.1*length(logf))),logs(round(0.8*length(logs))),alpha, ...
        'FontSize',12,'FontName','Helvetica');
    axis equal
    
    % manipulated noise spectra--------------------------------------------
    subplot(2,2,4);
    [f,s]=eo_polaraverage(abs(fft2(double(m_data))));
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
    alpha=sprintf('alpha=%2.2f',manipulated_im_slope);
    text(logf(round(0.1*length(logf))),logs(round(0.8*length(logs))),alpha,...
        'FontSize',12,'FontName','Helvetica');
    axis equal

end

end