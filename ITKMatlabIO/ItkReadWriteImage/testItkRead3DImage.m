% Test itkRead3DImage lib
% Author: Hui Xie
% date: Sep 20th, 2016

% 3D file: 'tenGlyphBqdAbcUv.nrrd'ls 'fmob-c4h.nrrd'

clear all;
imageArray = itkRead3DImage('fmob-c4h.nrrd');
%imageArray = itkRead3DImage('tenGlyphBqdAbcUv.nrrd');
imageArray = imageArray*100;

imageSize = size(imageArray)
for slice = 1:imageSize(3)
    image(imageArray(:,:,slice));
    fprintf('Slice: %d, ',slice);
    input('press return to next:');
end
close