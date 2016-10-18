%% Function to test the reading/writing of NRRD files
%% The purpose of this file is to test if nrrd files can
%% be read/written without loss of information
function [ errorStatus, outDWI, refDWI ] = testCS( )

    % Add main MATLAB_SCRIPTS code tree
    addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'))

    testFn = fullfile(fileparts(pwd),'TestSuite','Synthetic-Test-data_02.nhdr');
    testMaskFn = fullfile(fileparts(pwd),'TestSuite','mask.nhdr');
    refFn = fullfile(fileparts(pwd),'TestSuite','Synthetic-Ground-truth_02.nhdr');
    % Check that the file exists
    assert(exist(testFn, 'file') == 2, 'File does not exist');
    assert(exist(testMaskFn, 'file') == 2, 'File does not exist');
    assert(exist(refFn, 'file') == 2, 'File does not exist');

    outFn = fullfile(tempdir(),'csout.nhdr');

    tic
    %cache_result_fn = fullfile(tempdir(),'cache_03.mat');
    %if 1 ||  ~exist(cache_result_fn,'file')
        run_cs(testFn,testMaskFn,outFn, 0.01);
    %    save( cache_result_fn )
    %else
    %    load( cache_result_fn )
    %end
    toc

    outDWI = itkLoadWithMetadata(outFn) %nrrdLoadWithMetadata(outFn)
    refDWI = nrrdLoadWithMetadata(refFn) %nrrdLoadWithMetadata(refFn)

    errorStatus = compareDWIdata( refDWI, outDWI);
    if errorStatus < 0.05
       sprintf('SUCCESS: %s file matched',refFn)
    else
       sprintf('FAILED: %s file matched',refFn)
    end

    if comparevars(refDWI,outDWI,0.1)
        sprintf('SUCCESS: %s file matched',refFn)
    else
       sprintf('FAILED: %s file matched',refFn)
    end
end
