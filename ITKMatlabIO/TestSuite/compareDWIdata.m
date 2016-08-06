
function [ errorValue ] = compareDWIdata( dwi1, dwi2 )
    errorValue = 0.0;
    errorValue = errorValue + compareImages( dwi1, dwi2 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'bvalue',0.0, 0.0 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'gradientdirections',0.0001, 0.0001 );
    errorValue = errorValue + compareFloatArrays(dwi1, dwi2, 'measurementframe', 0.0, 0.0 );
end

function [ errorValue ]  = compareImages ( im1, im2 )
    errorValue = 0.0;
    errorValue = errorValue + abs(double(im1.space)-double(im2.space));
    errorValue = errorValue + compareFloatArrays(im1, im2, 'data', 0.015, 0.005 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'spacedirections', 0.0, 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'spaceorigin', 0.0, 0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'centerings', 0.0,  0.0 );
    errorValue = errorValue + compareFloatArrays(im1, im2, 'kinds',0.0, 0.0 );  % Check the enumerations are same

    errorValue = errorValue + compareStringArrays( im1, im2, 'spaceunits');
    errorValue = errorValue + compareStringArrays( im1, im2, 'spacedefinition');
end

function [ error ] = compareFloatArrays( fa, sa, fieldName, elementTolerance, cummulativeTolerance )
  %% Compare first and second array to determine how different they are
  if isfield(fa, fieldName) && isfield(sa, fieldName)
      array1 = double(fa.(fieldName));
      array2 = double(sa.(fieldName));
      denom = norm(array2(:))^2;
      if ~( isfinite(denom) && denom > eps )
          denom = 1;
      end
      error = norm(array1(:)-array2(:))^2/denom;
      if size(size(array1),2) == 4 % is data array
          size_voxels=size(array1,1)*size(array1,2)*size(array1,3);
          size_gradients=size(array1,4);
          array1=reshape(array1,[size_voxels,size_gradients]);
          array2=reshape(array2,[size_voxels,size_gradients]);
          errorPerElement=zeros(size_voxels,1);
          for i=1:size_voxels
            errorPerElement(i)=norm(array1(i,:) - array2(i,:))^2/norm(array1(i,:))^2;
          end
          errorElement = max(abs(errorPerElement));
          hist(errorPerElement);
      else
          errorElement = max(abs(array1(:) - array2(:)));
      end
  else
      sprintf('ERROR: MISSING FIELD: %s',fieldName);
      error = 10000000;
  end
  if ( ~isfinite(error) ) ||  (error > cummulativeTolerance)
      sprintf('ERROR: Differences found for: %s, (%f > %f)', fieldName, error, cummulativeTolerance)
      error = 10000000;
  end
  if errorElement > elementTolerance
      sprintf('ERROR: Differences found for element error: %s, (%f > %f)', fieldName, errorElement, elementTolerance)
      error = 10000000;
  end
end

function [ error ] = compareStringArrays( fa, sa, fieldName )
  %% Compare first and second array to determine how different they are
  if isfield(fa, fieldName) && isfield(sa, fieldName)
      sarray1 = fa.(fieldName);
      sarray2 = sa.(fieldName);
      error = sum(1.-strcmp(sarray1,sarray2));
  else
      sprintf('ERROR: MISSING FIELD: %s',fieldName);
      error = 10000000;
  end
  if error > 0
      sprintf('ERROR: Differences found for: %s', fieldName)
  end
end
