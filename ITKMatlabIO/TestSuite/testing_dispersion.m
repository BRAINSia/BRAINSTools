function [error_r1, error_r3, error_c1, error_c3] = testing_dispersion()
%testDispersion A single function for testing the code
% for computing dispersion

  base_path_slow_vtk_reader = fullfile(fileparts(pwd),'Medical_Image_Processing_Toolkit');

  addpath(fullfile(base_path_slow_vtk_reader,'class_image'))
  addpath(fullfile(base_path_slow_vtk_reader,'class_mesh'))
  addpath(fullfile(base_path_slow_vtk_reader,'IO'))

  addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'))

  % Determine if VTK or Mat files should be used
  useVTK = 0;  % <---- INFO: Peter This needs to be changed to a 1 and this test pass.
  useVTK = 1;
  for useVTK=1:1
      if useVTK
          testVTKRaysFn = fullfile(fileparts(pwd),'TestSuite','disp_test_rays.vtp');
          testVTKCirclesFn = fullfile(fileparts(pwd),'TestSuite','disp_test_circles.vtp');
          % Check that the file exists
          assert(exist(testVTKRaysFn, 'file') == 2, 'File does not exist');
          assert(exist(testVTKCirclesFn, 'file') == 2, 'File does not exist');
          RayStruct = vtkLoadPolyData(testVTKRaysFn);
          CircleStruct = vtkLoadPolyData(testVTKCirclesFn);
          rays = RayStruct.fibers;
          circles = CircleStruct.fibers;

          r1VTKFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_rays.vtp');
          r3VTKFn = fullfile(fileparts(pwd),'TestSuite','disp_scale3_numdirs10_1_1_rays.vtp');
          c1VTKFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_circles.vtp');
          c3VTKFn = fullfile(fileparts(pwd),'TestSuite','disp_scale3_numdirs10_1_1_circles.vtp');

          known_r1_struct = vtkLoadPolyData(r1VTKFn);
          known_r3_struct = vtkLoadPolyData(r3VTKFn);
          known_c1_struct = vtkLoadPolyData(c1VTKFn);
          known_c3_struct = vtkLoadPolyData(c3VTKFn);

          known_r1_struct.DDF = known_r1_struct.DDF';
          known_r3_struct.DDF = known_r3_struct.DDF';
          known_c1_struct.DDF = known_c1_struct.DDF';
          known_c3_struct.DDF = known_c3_struct.DDF';

          rmfield(known_r1_struct,'values')
          rmfield(known_r3_struct,'values')
          rmfield(known_c1_struct,'values')
          rmfield(known_c3_struct,'values')

          dispersionFeildName = 'DDF';
      else
          testMatRaysFn = fullfile(fileparts(pwd),'TestSuite','disp_test_rays.mat');
          testMatCirclesFn = fullfile(fileparts(pwd),'TestSuite','disp_test_circles.mat');
          % Check that the file exists
          assert(exist(testMatRaysFn, 'file') == 2, 'File does not exist');
          assert(exist(testMatCirclesFn, 'file') == 2, 'File does not exist');
          RayStruct = load(testMatRaysFn);
          CircleStruct = load(testMatCirclesFn);
          rays = RayStruct.disp_test_rays;
          circles = CircleStruct.disp_test_circles;

          r1MatFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_rays.mat');
          r3MatFn = fullfile(fileparts(pwd),'TestSuite','disp_scale3_numdirs10_1_1_rays.mat');
          c1MatFn = fullfile(fileparts(pwd),'TestSuite','disp_scale1_numdirs10_1_1_circles.mat');
          c3MatFn = fullfile(fileparts(pwd),'TestSuite','disp_scale3_numdirs10_1_1_circles.mat');

          known_r1_struct = load(r1MatFn,'DDF');
          known_r3_struct = load(r3MatFn,'DDF');
          known_c1_struct = load(c1MatFn,'DDF');
          known_c3_struct = load(c3MatFn,'DDF');
          dispersionFeildName = 'DDF';
      end

      NUMDIRS=10;
      testd_ddf_r1 = compute_dispersion(rays, 1, NUMDIRS, 1, 1);
      error_r1 = computeDispersionError ( testd_ddf_r1, known_r1_struct.(dispersionFeildName) )


      testd_ddf_r3 = compute_dispersion(rays, 3, NUMDIRS, 1, 1);
      error_r3 = computeDispersionError ( testd_ddf_r3, known_r3_struct.(dispersionFeildName) )

      testd_ddf_c1 = compute_dispersion(circles, 1, NUMDIRS, 1, 1);
      error_c1 = computeDispersionError ( testd_ddf_c1, known_c1_struct.(dispersionFeildName) )

      testd_ddf_c3 = compute_dispersion(circles, 3, NUMDIRS, 1, 1);
      error_c3 = computeDispersionError ( testd_ddf_c3, known_c3_struct.(dispersionFeildName) )
  end

  testWritingDecorated = 0;
  if testWritingDecorated
     known_r1_struct.(dispersionFeildName) = testd_ddf_r1';
     known_r3_struct.(dispersionFeildName) = testd_ddf_r3';
     known_c1_struct.(dispersionFeildName) = testd_ddf_c1';
     known_r3_struct.(dispersionFeildName) = testd_ddf_c3';

     vtkSavePolyData('test_r1_out.vtk',known_r1_struct);
     vtkSavePolyData('test_r3_out.vtk',known_r3_struct);
     vtkSavePolyData('test_c1_out.vtk',known_c1_struct);
     vtkSavePolyData('test_r3_out.vtk',known_c3_struct);
  end


end

function [ error ] = computeDispersionError( testd,refd )
  error_vector = abs(testd - refd);
  error = sum(error_vector);
end

function [] = writeMyDecoratedVTK(vtkFileName, vtkStructurePoint, dispersionDecorators)
  % writeMyDecoratedVTK This will write a Slicer compatible vtk file
  % that can be used to display or compute dispersion statistics across
  % all subjects.

  %  I don't know what is needed here to merge the dispersion info
  decoratedVTK = DecorateWithScalars( vtkStructurePoint, dispersionDecorators );
  write_vtkMesh(vtkFileName,decoratedVTK);
end

function [ decoratedVTK ] = DecorateWithScalars( vtkStructurePoint, dispersionDecorators )
  % DecorateWithScalars Given a vtkPolyData, and correctly sized scalars,
  % decorate the fiber tracts with the scalars.
  decoratedVTK = vtkStructurePoint;
end
