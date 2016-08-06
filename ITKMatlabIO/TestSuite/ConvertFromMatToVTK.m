
addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'))


PointDataRays=load('disp_test_rays.mat')
refRay=vtkLoadPolyData('disp_test_rays.vtk')
refRay.fibers=PointDataRays.disp_test_rays
vtkSavePolyData('disp_test_rays.vtp',refRay)


newRay=vtkLoadPolyData('disp_scale1_numdirs10_1_1_rays.vtk')
newRay.fibers=refRay.fibers
disp=load('disp_scale1_numdirs10_1_1_rays.mat')
[newRay.DDF]=disp.DDF'
rmfield(newRay,'values')
vtkSavePolyData('disp_scale1_numdirs10_1_1_rays.vtp',newRay)

newRay=vtkLoadPolyData('disp_scale3_numdirs10_1_1_rays.vtk')
newRay.fibers=refRay.fibers
disp=load('disp_scale3_numdirs10_1_1_rays.mat')
[newRay.DDF]=disp.DDF'
rmfield(newRay,'values')
vtkSavePolyData('disp_scale3_numdirs10_1_1_rays.vtp',newRay)

PointDataCircles=load('disp_test_circles.mat')
refCircle=vtkLoadPolyData('disp_test_circles.vtk')
refCircle.fibers=PointDataCircles.disp_test_circles
vtkSavePolyData('disp_test_circles.vtp',refCircle)

newCircle=vtkLoadPolyData('disp_scale1_numdirs10_1_1_circles.vtk')
newCircle.fibers=refCircle.fibers
disp=load('disp_scale1_numdirs10_1_1_circles.mat')
[newCircle.DDF]=disp.DDF'
rmfield(newCircle,'values')
vtkSavePolyData('disp_scale1_numdirs10_1_1_circles.vtp',newCircle)

newCircle=vtkLoadPolyData('disp_scale3_numdirs10_1_1_circles.vtk')
newCircle.fibers=refCircle.fibers
[newCircle.DDF]=disp.DDF'
rmfield(newCircle,'values')
vtkSavePolyData('disp_scale3_numdirs10_1_1_circles.vtp',newCircle)

