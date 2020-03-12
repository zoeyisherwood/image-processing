function [output_im] = change_rms_contrast(input,desired_contrast,plotfigs)
%[output_im] = change_rms_contrast(input,desired_contrast,plotfigs)
%
% the bulk of this code works off the lumMatch.m function from the SHINE
% toolbox (Verena Willenbockel -
% http://www.mapageweb.umontreal.ca/gosselif/SHINE/)
% 
%   input - your input image loaded in using imread
%   desired_contrast - input the output contrast you want (range: 0 to 1)
%        note that the higher the contrast, the higher the likelihood of
%        clipping.
%   plotfigs - set to 1 if you want output figures of resultant image and
%               histogram
% 
% program will change the rms contrast of your input image based on your
% desired contrast level (desired_contrast)
% 
% the way this code modifies contrast is by keeping mean luminance the same
% and then changing the standard deviation of pixels to yield an rms
% contrast level of your choosing.
% 
%  equation:
% 
%       rmscontrast = std2(input(:))/mean(input(:)
% 
%       desiredrmscontrast = X / mean(input(:))
%           X = desiredrmscontrast*meaninput(:)
%
%       **set std of input to X to get desired rms contrast.
% 
% current limitations of script:
%   -converts rgb input into grayscale. plan to add compatibility
%   for color images...
%   -data only output in uint8 format, which can cause clipping artefacts
%   for steep slopes. plan to add higher bitrate options (2^10, 2^14) which
%   can be used on high bitrate monitors like Display++...
%
% dependent scripts:
%   you need to have:
%       -lumMatch.m (http://www.mapageweb.umontreal.ca/gosselif/SHINE/)
% 
% log:
% 20200310: Initialised. zoeyisherwood.
% contact: zoey.isherwood@gmail.com

%% start processing...

% check if rgb-------------------------------------------------------------

if ndims(input) == 3
    disp(['your input is colored. your input will be' ...
        ' converted to greyscale...'])
    input = input(:,:,1:3); %just in case there's an alpha channel in there.
    input = rgb2gray(input);
    
end

% report rms contrast of input---------------------------------------------

orig_rms = std2(input(:))/mean(input(:));

disp(['input rms contrast: ' num2str(orig_rms)]);
disp(['desired rms contrast: ' num2str(desired_contrast)]);

% calculate what luminance stdev we need to get the desired contrast we
% want---------------------------------------------------------------------

stdev = desired_contrast*mean(input(:));

% use lumMatch.m to manipulate RMS contrast--------------------------------

output_im = lumMatch({input},[],[mean(input(:)) stdev]);
output_im = output_im{1};

% measure rms of manipulated image-----------------------------------------

output_rms = std2(output_im(:))/mean(output_im(:));


%plot figs if indicated----------------------------------------------------

if plotfigs == 1
    
    figure;
    
    % plot original image--------------------------------------------------
    subplot(2,2,1)
    imshow(input)
    title('input im')
    
    % plot output image----------------------------------------------------
    subplot(2,2,2)
    imshow(output_im)
    title('output im')
    
    % plot histogram of orig image-----------------------------------------
    subplot(2,2,3)
    hist(single(input(:)))
    xlim([0 255])
    title(['rms contrast: ' num2str(orig_rms) ' mean: ' num2str(mean(input(:)))]);
    
    % plot histogram of manipulated image---------------------------------
    subplot(2,2,4)
    hist(single(output_im(:)))
    xlim([0 255])
    title(['rms contrast: ' num2str(output_rms) ' mean: ' num2str(mean(output_im(:)))]);
    
  
end


end