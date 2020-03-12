function [outputim,maskmap] = makeimagefactortwo(input)
% pads your image with zeros if it's not square and ensures its dims are a
% factor of 2.
% eg, a 390x512 image would become 512x512
%
% this code will also spit out a mask map so the image can be easily
% cropped after analysis.
%
% can use maskmap to crop image back to original size with following lines;
%
%   [row_start cols_start] = find(maskmap, 1, 'first')
%   [row_end cols_end] = find(maskmap, 1, 'last')
%   finalim = outputim(row_start:row_end,cols_start:cols_end);
%
% Log:
% 20200311. Initialised. zoeyisherwood.
% Contact: zoey.isherwood@gmail.com

% check dims---------------------------------------------------------------

if ndims(input) == 3
    disp(['your input is colored. your input will be' ...
        ' converted to greyscale...'])
    input = input(:,:,1:3); %just in case there's an alpha channel in there.
    input = rgb2gray(input);
    
end

% get sizes----------------------------------------------------------------

a = size(input,1);
b = size(input,2);

% check if dimensions are a factor of 2------------------------------------

factortwovalue = max(2 .^ ceil(log2(size(input))));


if a == factortwovalue && b == factortwovalue %end function if the image is already square
    
    disp('your input is already square. and the dims are a factor of two')
    outputim = input;
    maskmap = ones(a,b);
    
elseif a ~= factortwovalue || b ~= factortwovalue
    
    disp('does not have dims that are a factor of two. converting...');
    
    maskmap = ones(a,b);
    
    %determine padding amount----------------------------------------------
    
    %firstly, which dim is larger?
    
    padsize = max(factortwovalue,factortwovalue);
    
    % ... and is this dim an even number? if not, add 1.
    
    % pad depending on dimension-------------------------------------------
    
    outputim = zeros(padsize,padsize);
    
    start_row = ceil((padsize-a)/2 + 1);
    start_col = ceil((padsize-b)/2 + 1);
    
    outputim(start_row:start_row+a-1,start_col:start_col+b-1) = input;
    
    maskmap = zeros(padsize,padsize);
    maskmap(start_row:start_row+a-1,start_col:start_col+b-1) = ones(a,b);
    
    disp('done')
    
end

end
