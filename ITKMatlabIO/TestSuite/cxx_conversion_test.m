function [error_r1] = cxx_conversion_test()
% cxx_conversion_test A single function for testing the code
% for computing dispersion

% The purpose of this function is to facilitate conversion to c++
% The matlab features will be replaced by eigen where possible.
% http://eigen.tuxfamily.org/dox-devel/AsciiQuickReference.txt

  % Like a #include for c++
  addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'))
  % Set up file paths to test data
  testVTKRaysFn = fullfile(fileparts(pwd),'TestSuite','disp_test_rays.vtk');
  % Check that the file exists
  assert(exist(testVTKRaysFn, 'file') == 2, 'File does not exist');
  % Load up some test data, and format for testing  <-- This part will
  % eventually be part of the c++ testing suite as well
  r1VTKFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_rays.vtk');
  known_r1_struct = vtkLoadPolyData(r1VTKFn);
  known_r1_struct.DDF = known_r1_struct.DDF';
  rmfield(known_r1_struct,'values')
  dispersionFeildName = 'DDF';

  %%%%%%%%%%%%%
  % HERE HERE HERE -- Start paying close attention here the
  % following section needs to converted to c++
  % fiberdispersion --inputFiberBundle testVTKRaysFn \
  %                 --outputFiberBundle decoratedOutput.vtk \
  %                 --scale <defaults to 1> \
  %                 --numDirs <defaults to 10> \
  %                 --tractSubsampling <defaults to 1> \
  %                 --fiberPointSubsampling <defaults to 1> \
  % // DO NOT IMPLEMENT --debugOutputFile <-- Ignore this
  RayStruct = vtkLoadPolyData(testVTKRaysFn);
  rays = RayStruct.fibers;
  %There are 22 fibers in this file, each fiber having exactly 11 points
  %RayStruct =
  %   fibers: {1x22 cell}   <- {1x numFibers cell}
  %RayStruct.fibers
  %  [3x11 double] ...  [3x11 double]  {3xnumPointsinFiber} where the 3
  %  points represent x,y,z physical space.
  % NOTE: For this particular case, there are 11 points per fiber for 22
  % fibers = 242 total points.


  NUMDIRS=10;
  testd_ddf_r1 = compute_dispersion(rays, 1, NUMDIRS, 1, 1);
  %At this point, testd_ddf = 1x242 double corresponding to the points in
  %all the fibers.  The first 11 values belong to fiber1, the last 11
  %valuse belong to fiber22.
  [RayStruct.DDF]=testd_ddf_r1'  %NOTE:  Values here need transposing for writing purposes.
  vtkSavePolyData('decoratedOutput.vtk',RayStruct)

  error_r1 = computeDispersionError ( testd_ddf_r1, known_r1_struct.(dispersionFeildName) )

end

function [ error ] = computeDispersionError( testd,refd )
  error_vector = abs(testd - refd);
  error = sum(error_vector);
end

