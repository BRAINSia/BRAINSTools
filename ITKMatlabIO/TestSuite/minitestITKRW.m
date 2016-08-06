%% TEST DWI load/save from ITKv4 for files previously written by nrrdReadWrite mechanisms

    %% NOTE: It appears that when a valid DWI file is provided,
    %%       That both the itkLoadWithMetadata and nrrdLoadWithMetadata functions perform
    %%       similarly
    inFn = fullfile(fileparts(pwd),'TestSuite','nrrd_writer_out.nhdr');
    inFn
    inDWIfromITK=itkLoadWithMetadata(inFn);
    inDWIfromNRRD=nrrdLoadWithMetadata(inFn);
    inDWIfromITK
    inDWIfromNRRD


    outFn = fullfile(tempdir(),'FromITKDWI.nrrd');%nrrd, not FromITKDWI.nhdr
    outFn

    itkSaveWithMetadata(outFn,inDWIfromITK); %% NOTE THIS  APPEARS TO PRODUCED AN INVALID DWI FILE
    %% cat outFn
%% NRRD0005
%%# Complete NRRD file format specification at:
%%# http://teem.sourceforge.net/nrrd/format.html
%%type: double
%%dimension: 4
%%space dimension: 4
%%sizes: 70 20 2 82
%% ^^^^^^^^^^^^^^^^^^^^^^
%% NOTE:  This is a 4D image rather than a 3D + gradients image!
    %%% KENT The following line is causing a test failure.
    inDWIAfterSaveITK=itkLoadWithMetadata(outFn);
    inDWIAfterSaveITK1=nrrdLoadWithMetadata(outFn);%testline
    inDWIAfterSaveITK
    inDWIAfterSaveITK1
