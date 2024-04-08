% Goal:build GroudTruth for itkLoadWithMetadata and itkSaveWithMetadata
% Function: Read all nhdr files in specific dir, and save their mat file in same name
% Author: Hui Xie
% Date: Oct 11th, 2016
% notes:
%    1 this program must locate outside GroundTruth dir;
%    2 after new program, please DON'T run this proram, which will
%      overwrite exising mat files;

disp('=============');
disp('Warning: After new program, please don''t run this program, which will overwrite the existing mat files.');
choice = input('Are you sure to run this program ? (yes/no):','s');
if 0 == strcmpi(choice,'yes')
    disp('=====Program exit=====')
    return;
end

clear all;
dirName = './GroundTruth/';
fileList = dir(dirName);
fileListSize = length(fileList);


disp('Load files:');
num = 0;
for i = 1:fileListSize
    if length(fileList(i).name) < 4
        continue; % ignore the . and .. dir files
    end

    [~,name,ext] = fileparts(fileList(i).name);
    if 0 == strcmp(ext,'.nhdr')
        continue  % ignore non-hhdr files
    end

    pathFilename = strcat(dirName,fileList(i).name);
    disp(pathFilename);
    loadStruct = itkLoadWithMetadata(pathFilename);
    matFilename = strcat(dirName,name,'.mat');
    delete(matFilename); % if file exist, delete it
    save(matFilename,'loadStruct');
    num =num+ 1;
end
fprintf('\nTotal %d nhdr files loaded and saved into %s dir in mat format\n', num, dirName);
disp('======Program End=======')
