function [ errorStatus, vtkStruct ] = testVTKReading( )
%testVTKReading Summary of this function goes here
%   Test reading of vtk files

  % Add main MATLAB_SCRIPTS code tree
  addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'))

  base_path_slow_vtk_reader = fullfile(fileparts(pwd),'Medical_Image_Processing_Toolkit');

  addpath(fullfile(base_path_slow_vtk_reader,'class_image'))
  addpath(fullfile(base_path_slow_vtk_reader,'class_mesh'))
  addpath(fullfile(base_path_slow_vtk_reader,'IO'))

  testFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_circles.vtk');
  % Check that the file exists
  assert(exist(testFn, 'file') == 2, 'File does not exist');

  vtkStruct = read_vtkMesh(testFn);

  % KENT TODO:  Write a new binary vtk mex reader that emulates the
  % behavior of the slow "read_vtkMesh" function.
  % As an added bonus, the binary reader should be able to
  % vtkStruct = vtkReadPolyData();

  % KENT TODO: Write a vtk polydata binary mex writer
  % Necesssary requirement:  the binary vtk readers/writers should
  % be able to read a "decorated" UKF generated .vtp file and write
  % it back out with loss of information.

  errorStatus = 1000; %% FAILED code
end

