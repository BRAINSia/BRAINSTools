% Goal:Compare GroudTruth with new itkLoadWithMetadata and itkSaveWithMetadata
% Function:
%   1  Using new itkLoad* and itkSave* programs, read all nhdr files in specific dir;
%   2  and compare them with exist mat files;
%   3  if they are inconsistent, give out information;
%   4  this program does not output any mat files in order to avoid to
%      overwriting existing ground truth mat files;
% notes: this program must locate outside GroundTruth dir;
% Author: Hui Xie
% Date: Oct 11th, 2016

clear all;
dirName = './GroundTruth/';
fileList = dir(dirName);
fileListSize = length(fileList);
disp('=========================');
disp('Compare GroudTruth with new itkLoadWithMetadata and itkSaveWithMetadata');
disp('Load files to check whether it is consistent with existing mat file:');
disp('  ');
totalNum = 0;
diffNum = 0;
for i = 1:fileListSize
    if length(fileList(i).name) < 4
        continue; % ignore the . and .. dir files
    end

    [~,name,ext] = fileparts(fileList(i).name);
    if 0 == strcmp(ext,'.nhdr')
        continue  % ignore non-hhdr files
    end

    pathFilename = strcat(dirName,fileList(i).name);
    loadStruct = itkLoadWithMetadata(pathFilename);

    % ===========For Test start
    %loadStruct.data(200) = 900;
    %loadStruct = rmfield(loadStruct,'kinds');
    %loadStruct.spacedefinition = 'Love this game';
    %loadStruct.spaceunits{2} = 'cm';
    % ===========For Test end


    matFilename = strcat(dirName,name,'.mat');
    existingStruct = importdata(matFilename);


    totalNum =totalNum+ 1;
    if 1 == isequal(loadStruct, existingStruct)
        fprintf('%s : is consistent.\n', pathFilename);
    else
        disp(pathFilename);
        % check fieldnames
        loadFieldNames = fieldnames(loadStruct);
        existFieldNames = fieldnames(existingStruct);
        if 0 == isequal(loadFieldNames, existFieldNames)
            disp('========Inconsistent fieldname:');
            fprintf('loadStruct has %d fields:\n', length(loadFieldNames));
            disp(loadFieldNames);
            fprintf('existingStruct has %d fields:\n', length(existFieldNames));
            disp(existFieldNames);
            return;
        end

        % check field element
        numFields = length(loadFieldNames);
        for i = 1:numFields
            fieldName = loadFieldNames{i};
            if 0 == isequal(loadStruct.(fieldName), existingStruct.(fieldName))
                disp(strcat('==========Inconsistent in field ####', fieldName, '####'));
                nL = numel(loadStruct.(fieldName));
                nE = numel(loadStruct.(fieldName));
                if (nL ~= nE)
                    fprintf('=========different number of elements: LoadStruct has %d, while exisingStruct has %d elements;\n', nL, nE);
                    return
                else
                    if isnumeric(loadStruct.(fieldName))
                        for j=1:nL
                            if loadStruct.(fieldName)(j) ~= existingStruct.(fieldName)(j)
                                fprintf('==========in element %d, LoadStruct has %d, while ExisingStruct has %d;\n', j,loadStruct.(fieldName)(j), existingStruct.(fieldName)(j));
                                return;
                            end
                        end
                    elseif iscell(loadStruct.(fieldName))
                        for j=1:nL
                            if 0 == isequal(loadStruct.(fieldName)(j),existingStruct.(fieldName)(j))
                                fprintf('==========in element %d, LoadStruct has ''%s'', while ExisingStruct has ''%s'';\n', j,loadStruct.(fieldName){j}, existingStruct.(fieldName){j});
                                return;
                            end
                        end
                    else
                        if 0 == strcmp(loadStruct.(fieldName),existingStruct.(fieldName))
                            fprintf('==========in element %d, LoadStruct has ''%s'', while ExisingStruct has ''%s'';\n', j,loadStruct.(fieldName), existingStruct.(fieldName));
                            return;
                        end
                    end

                end

            end
        end

    end


end
fprintf('\nTotally verify %d nhdr files with existing mat files. They are consistent. !!!Good!!!.\n', totalNum);
disp('======Program End=======')