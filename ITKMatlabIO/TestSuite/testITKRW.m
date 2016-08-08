%% Function to test the reading/writing of NRRD files
%% The purpose of this file is to test if nrrd files can
%% be read/written without loss of information
function [ ] = testITKRW( )
    % Add main MATLAB_SCRIPTS code tree
    addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'));

    nrrd4DFloat = fullfile(fileparts(pwd),'TestSuite','nrrd_writer_out.nhdr');
    nrrd3DShort = fullfile(fileparts(pwd),'TestSuite','mask.nhdr');

    % Check that the file exists
    assert(exist(nrrd4DFloat, 'file') == 2, 'File does not exist');
    assert(exist(nrrd3DShort, 'file') == 2, 'File does not exist');

    %======================================================================
    %% TEST DWI load/save from ITKv4 for files previously written by nrrdReadWrite mechanisms (4D Float NRRD)

    itkLoad4DFloat=itkLoadWithMetadata(nrrd4DFloat);
    nrrdLoad4DFloat=nrrdLoadWithMetadata(nrrd4DFloat);

%    nrrdLoad4DFloat.spaceorigin = nrrdLoad4DFloat.spaceorigin'
%    nrrdLoad4DFloat.spaceorigin = [1,2,3]
%    nrrdLoad4DFloat;
%    ITKimage = MatlabITKImageObject(itkLoad4DFloat);
%    ITKimage.PrintSelf()
%    ITKimage.ToMatlabStruct()


    if comparevars(nrrdLoad4DFloat,itkLoad4DFloat);
       sprintf('SUCCESS: nrrd4DFloat loading with nrrd and itk')
    else
       sprintf('FAILED: nrrd4DFloat loading with nrrd and itk')
    end

    %======================================================================
    %% TEST DWI load/save from ITKv4 for files previously written by nrrdReadWrite mechanisms (3D Short NRRD)

    itkLoad3DShort=itkLoadWithMetadata(nrrd3DShort);
    nrrdLoad3DShort=nrrdLoadWithMetadata(nrrd3DShort);

    if comparevars(itkLoad3DShort,nrrdLoad3DShort);
    sprintf('SUCCESS: nrrd3DShort loading with nrrd and itk')
    else
    sprintf('FAILED: nrrd3DShort loading with nrrd and itk')
    end

    %======================================================================
    %% TEST if saving and loading from ITK is consistent. (4D Float NRRD)
    itkSave4DFloat = fullfile(tempdir(),'itkSave4DFloat.nrrd');
    itkSaveWithMetadata(itkSave4DFloat,itkLoad4DFloat);
    itkLoad4DFloatItkSave4DFloat=itkLoadWithMetadata(itkSave4DFloat);


    if comparevars(itkLoad4DFloat,itkLoad4DFloatItkSave4DFloat);
       sprintf('SUCCESS: itkLoad4DFloat,itkLoad4DFloatItkSave4DFloat')
    else
       sprintf('FAILED: itkLoad4DFloat,itkLoad4DFloatItkSave4DFloat')
    end

    %======================================================================
    %% TEST if saving and loading from nrrd is consistent. (4D float NRRD)
    nrrdSave4DFloat = fullfile(tempdir(),'nrrdSave4DFloat.nrrd');
    nrrdSaveWithMetadata(nrrdSave4DFloat,nrrdLoad4DFloat);
    nrrdLoad4DFloatNrrdSave4DFloat=nrrdLoadWithMetadata(nrrdSave4DFloat);


    if comparevars(nrrdLoad4DFloat,nrrdLoad4DFloatNrrdSave4DFloat);
    sprintf('SUCCESS: nrrdLoad4DFloat,nrrdLoad4DFloatNrrdSave4DFloat')
    else
    sprintf('FAILED: nrrdLoad4DFloat,nrrdLoad4DFloatNrrdSave4DFloat')
    end

    %======================================================================
    %% TEST if final RW from itk and nrrd are comparable (4D float nrrd)
    errorStatus = comparevars(itkLoad4DFloatItkSave4DFloat,nrrdLoad4DFloatNrrdSave4DFloat);
    if errorStatus;
    sprintf('SUCCESS: itkLoad4DFloatItkSave4DFloat,nrrdLoad4DFloatNrrdSave4DFloat')
    else
    sprintf('FAILED: itkLoad4DFloatItkSave4DFloat,nrrdLoad4DFloatNrrdSave4DFloat')
    end

    %======================================================================
    %% TEST if saving and loading from ITK is consistent. (3D short NRRD)
    itkSave3DShort = fullfile(tempdir(),'itkSave3DShort.nrrd');
    itkSaveWithMetadata(itkSave3DShort,itkLoad3DShort);
    itkLoad3DShortItkSave3DShort=itkLoadWithMetadata(itkSave3DShort);

    errorStatus = comparevars(itkLoad3DShort,itkLoad3DShortItkSave3DShort);
    if errorStatus;
    sprintf('SUCCESS: itkLoad3DShort,itkLoad3DShortItkSave3DShort')
    else
    sprintf('FAILED: itkLoad3DShort,itkLoad3DShortItkSave3DShort')
    end

    %======================================================================
    %% TEST if saving and loading from nrrd is consistent. (3D short NRRD)
    nrrdSave3DShort = fullfile(tempdir(),'nrrdSave3DShort.nrrd');
    nrrdSaveWithMetadata(nrrdSave3DShort,nrrdLoad3DShort);
    nrrdLoad3DShortNrrdSave3DShort=nrrdLoadWithMetadata(nrrdSave3DShort);

    errorStatus = comparevars(nrrdLoad3DShort,nrrdLoad3DShortNrrdSave3DShort);
    if errorStatus;
    sprintf('SUCCESS: nrrdLoad3DShort,nrrdLoad3DShortNrrdSave3DShort')
    else
    sprintf('FAILED: nrrdLoad3DShort,nrrdLoad3DShortNrrdSave3DShort')
    end

    %======================================================================
    %% TEST if final RW from itk and nrrd are comparable (3D short NRRD)
    errorStatus = comparevars(itkLoad3DShortItkSave3DShort,nrrdLoad3DShortNrrdSave3DShort);
    if errorStatus;
    sprintf('SUCCESS: itkLoad3DShortItkSave3DShort,nrrdLoad3DShortNrrdSave3DShort')
    else
    sprintf('FAILED: itkLoad3DShortItkSave3DShort,nrrdLoad3DShortNrrdSave3DShort')
    end


end

% function [] = dummyCommentOut()
%
%     %% cat outFn
% %% NRRD0005
% %%# Complete NRRD file format specification at:
% %%# http://teem.sourceforge.net/nrrd/format.html
% %%type: double
% %%dimension: 4
% %%space dimension: 4
% %%sizes: 70 20 2 82
% %% ^^^^^^^^^^^^^^^^^^^^^^
% %% NOTE:  This is a 4D image rather than a 3D + gradients image!
%     % KENT The following line is causing a test failure.
%     inDWIAfterSaveITK=itkLoadWithMetadata(outFn);
%     inDWIAfterSaveITK
%
%
%     %% Test writing in different formats
%     inDWI = itkLoadWithMetadata(testFn);
%
%     %% float
%     outDWI = inDWI
%     outDWI.data = single(inDWI.data);
%     outFn = fullfile(tempdir(),'test_single.nhdr');
%     outFn
%     itkSaveWithMetadata(outFn,outDWI);
%
%
%     %% double
%     outDWI = inDWI
%     outDWI.data = double(inDWI.data);
%     outFn = fullfile(tempdir(),'test_double.nhdr');
%     outFn
%     itkSaveWithMetadata(outFn,outDWI);
%
%    %% int32
%     outDWI = inDWI
%     outDWI.data = int32(inDWI.data);
%     outFn = fullfile(tempdir(),'test_int32.nhdr');
%     outFn
%     itkSaveWithMetadata(outFn,outDWI);
%
% end
%
