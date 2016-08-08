%% Function to test the reading/writing of NRRD files
%% The purpose of this file is to test if nrrd files can
%% be read/written without loss of information
function [ ] = TestITKOO( )
    % Add main MATLAB_SCRIPTS code tree
    addpath(fullfile(fileparts(pwd),'MATLAB_SCRIPTS'));

    nrrd4DFloat = fullfile(fileparts(pwd),'TestSuite','nrrd_writer_out.nhdr');
    HDnrrd = fullfile(fileparts(pwd),'TestSuite','example_hd_data.nrrd');
    HDnrrd = ('/Users/ericpahl/ResearchJohnson/CompressedSensing/TestSuite/example_hd_data.nrrd');

    nrrd3DShort = fullfile(fileparts(pwd),'TestSuite','mask.nhdr');

    % Check that the file exists
    assert(exist(nrrd4DFloat, 'file') == 2, 'File does not exist');
    assert(exist(nrrd3DShort, 'file') == 2, 'File does not exist');

    %======================================================================
    %% TEST DWI load/save from ITKv4 for files previously written by nrrdReadWrite mechanisms (4D Float NRRD)

    itkLoad4DFloat=itkLoadWithMetadata(nrrd4DFloat);
    nrrdLoad = nrrdLoadWithMetadata(HDnrrd);
    %TODO ERIC: eventually remove all the loading above and just have
    %manually declared structure

    obj= MatlabITKImageObject(nrrd4DFloat); % try NRRD Filename constructor
%     obj2 = MatlabITKImageObject(itkLoad4DFloat); % try NRRD Struct constructor
%     obj3 = MatlabITKImageObject(HDnrrd); %try NRRD struct constructor
%     naked = MatlabITKImageObject(); % try default constructor
     naked2 = MatlabITKImageObject(uint64(2)); %try 2D default constructor
%     copy = MatlabITKImageObject(obj);% try copy constructor

%     if obj==obj2 && obj2==copy && copy==obj3
%         disp( 'all objects are the same')
%     else
%         disp( 'objects are NOT the same')
%     end

%     obj.SetDirection([1 0 2;0 3 0;4 0 5]); %for testing inverse direction

    obj.PrintSelf();



%     naked.PrintSelf();
%
%     naked.SetOrigin([10;20;30]) ;
%
%     naked.PrintSelf();

%     foo(naked)%passing my simulated  reference
    %handles are copied by value as pointers
    %

%     ITKnrrd = FromNRRDStruct(itkLoad4DFloat);
%     ITKnrrd.PrintSelf();




%     SetOrigin(naked, [10;20;30]);
%
%       naked.PrintSelf();
%
%     disp( naked.GetOrigin() )
%
%     ITKimage = MatlabITKImageObject(itkLoad4DFloat);
%     ITKimage.PrintSelf()
%     ITKimage.ToMatlabStruct()
%
%
end


function foo(mo)
   mo.SetOrigin([15;25;35]);
end
