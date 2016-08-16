function [normArr,oldMin,scaleFactor] = NormalizeDataComponent(arr)
  % This function normalizes a 3D matrix between zero and one.
  newMax = single(1);
  newMin = single(0);
  oldMax = single(max(arr(:)));
  oldMin = single(min(arr(:)));
  scaleFactor = (newMax-newMin)/(oldMax-oldMin);
  normArr = (single(arr)-oldMin)*scaleFactor+newMin;
  normArr = single(normArr);
end
