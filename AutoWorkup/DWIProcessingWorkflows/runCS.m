function runCS(testFileName,maskFileName,outFileName)

    % Add main MATLAB_SCRIPTS code tree
    addpath('/Shared/sinapse/scratch/aghayoor/CompressedSensing/WorkInProgress/MATLAB_SCRIPTS')
    addpath('/Shared/sinapse/scratch/aghayoor/CompressedSensing/WorkInProgress/mexFiles')

    % Check that the file exists
    assert(exist(testFileName, 'file') == 2, 'File does not exist');
    assert(exist(maskFileName, 'file') == 2, 'File does not exist');

    tic
        run_cs(testFileName,maskFileName,outFileName, 0.01);
    toc

end
