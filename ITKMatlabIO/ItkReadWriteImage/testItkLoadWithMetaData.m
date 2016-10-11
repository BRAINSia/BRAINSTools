% Test loadWithMeatData file
% Author: Hui Xie
% date: Sep 20th, 2016

% 3D file: 'tenGlyphBqdAbcUv.nrrd'ls 'fmob-c4h.nrrd'

clear all;
imageStructure = itkLoadWithMetadata('fmob-c4h.nrrd')

imageArray = imageStructure.data*100;

imageSize = size(imageArray)
for slice = 1:imageSize(3)
    image(imageArray(:,:,slice));
    fprintf('Slice: %d, ',slice);
    input('press return to next:');
end
close