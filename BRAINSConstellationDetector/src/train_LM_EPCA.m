% Author: Wei Lu
% at Psychiatry Imaging Lab, University of Iowa Health Care, 2010
%
% This is the implementation of the two training phases for the proposed
% LME-EPCA algorithm (Section 3.3 in my MS thesis).
%
% Objective:
% Train a linear model that can evolutionary estimate new landmarks from
% already known landmarks.
%
% Input:
% baselandmarks     - Base landmarks in training datasets
% newLandmarks      - EPCA landmarks in training datasets
% numBaseLandmarks  - Number of base landmarks
% numNewLandmarks   - Number of EPCA landmarks
% numDatasets       - Number of training datasets
% dim               - Dimension of the landmarks
%
% Output: (variables stored in files)
% M                 - Optimal linear combination of principal components
%                     in each iteration
% s                 - The 3-by-1 mean of landmark vector space in each
%                     iteration
% search_radius     - Search radius of each EPCA landmark
% newLandmarksNames  - A name list of all EPCA landmarks


%% Dependencies
load_landmarks % load all landmarks and other input parameters


%% Init
W = cell(numNewLandmarks, 1);
s = zeros(dim, numNewLandmarks);
M = cell(numNewLandmarks, 1);
lmk_est = zeros(dim,numNewLandmarks,numDatasets);
err = zeros(dim,numNewLandmarks,numDatasets); % to find search_radius

% Initialize the landmark vector space
Xi = zeros(dim*(numBaseLandmarks-1), numDatasets);
for k = 2:numBaseLandmarks
    lmkVec = (baseLandmarks{k} - baseLandmarks{1})';
   Xi((k-2)*dim + 1 : (k-1)*dim, :) = lmkVec;
end

% Compute all principal components

% Test stuff
ratioPC1 = zeros(numNewLandmarks, 1);
ratioPC = zeros(numNewLandmarks, 1);

for i = 1:numNewLandmarks
    %% Training phase-1
    % Train principal components for landmark vector space associated with the
    % EPCA landmarks in each iteration.
    if i > 1
        lmkVec = (newLandmarks{i-1} - baseLandmarks{1})';
        Xi = [Xi; lmkVec];
    end

    % Remove Xi mean
    Xi_mean = zeros(dim, 1);
    mean_stacked = mean(Xi, 2);
    for d = 1:dim
        Xi_mean(d) = sum(mean_stacked(...
            d:dim:end-dim+d))/(numBaseLandmarks+i-2);
    end
    s(:, i) = Xi_mean;
    I_si = repmat(s(:, i), numBaseLandmarks+i-2, numDatasets);
    Xi_demeaned = Xi - I_si;
    [V, D] = eig(Xi_demeaned*Xi_demeaned');

    ratioPC1(i) = sum(sum(D(:, end))) / sum(D(:));

    % Number of PCs should be chosen so that the following condition is
    % met:
    % sum(sum(D(:, end-numPCs+1))) / sum(D(:)) > 99%;
    % in this case, numPCs >= 1 will meet the requirement
    % or we can use as must PCs as possible
    % tol argument of the rank func is set conservatively
    % Adjust the tol argument for different training datasets
    %tolXi = 100;
    %numPCs = rank(Xi_demeaned, tolXi);

    numPCs = 10;
    W{i} = V(:, end-numPCs+1 : end);

    %% Training phase-2
    % Train optimal linear relationship between already known landmark
    % vector space and the EPCA landmark to be estimated in each iteration.
    Zi = W{i}'*Xi_demeaned; % PCA mapped space
    Yi = newLandmarks{i} - baseLandmarks{1};
    Ci = (Zi*Zi')\Zi*Yi;
    M{i} = W{i}*Ci;

    %% Compute the estimation errors for training datasets
    Xi_t = zeros((numBaseLandmarks+i-2)*dim, 1);
    x1_t = baseLandmarks{1};
    I_si_t = repmat(s(:, i), numBaseLandmarks+i-2, 1);
    for j = 1:numDatasets
        for k = 2:numBaseLandmarks+i-1
            if k <= numBaseLandmarks
                xk_t = baseLandmarks{k};
            else
                xk_t = newLandmarks{k - numBaseLandmarks};
            end
            Xi_t((k-2)*dim+1:(k-1)*dim, 1) = xk_t(j, :)' - x1_t(j, :)';
        end
        Xi_demeaned_t = Xi_t - I_si_t;
        lmk_est(:, i, j) = x1_t(j, :)' + M{i}'*Xi_demeaned_t;
        lmk_groundTruth = newLandmarks{i};
        err(:, i, j) = lmk_est(:, i, j) - lmk_groundTruth(j, :)';
    end
end


%% Determine search radii
err_dist = zeros(numNewLandmarks,numDatasets); % Euclidean err distance
for j = 1:numNewLandmarks
    for k = 1:numDatasets
        err_dist(j,k) = sqrt(sum(err(:,j,k).^2));
    end
end
err_mean = mean(err_dist, 2);
err_std = sqrt(sum((err_dist - repmat(err_mean,1,numDatasets)).^2, ...
    2)/numDatasets);
err_max = max(err_dist,[],2);
err_min = min(err_dist,[],2);
err_median = median(err_dist,2);

% Aim for 99.7% coverage of anotomical variation (assuming normal distr.)
search_radius = err_mean + 3*err_std;
search_radius_max = 1.2*err_max;
search_radius_min = 1.6*ones(numNewLandmarks, 1); % if max < min, take min
search_radius(search_radius > search_radius_max) = ...
    search_radius_max(search_radius > search_radius_max);
search_radius(search_radius < search_radius_min) = ...
    search_radius_min(search_radius < search_radius_min);


%% Write model parameters to file
% txt version
fid = fopen('LME_EPCA.txt','w');
for i = 1:numNewLandmarks
    fprintf(fid, '%s\n', newLandmarksNames{i});
    for d = 1:dim
        fprintf(fid, '%6.32f ', s(d,i));
    end
    fprintf(fid, '\n');
    fprintf(fid, '%6.32f\n', search_radius(i));
    fprintf(fid, '%d\n', size(M{i},1));
    tmp = M{i};
    for d = 1:dim
        fprintf(fid, '%6.32f ', tmp(:,d));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
end
fclose(fid);

% mat version
fid = fopen('processingList.txt','w');
for i = 1:numNewLandmarks
    name = newLandmarksNames{i};
    fprintf(fid, '%s\n', name);
    fprintf(fid, '%6.32f\n', search_radius(i));
end
fclose(fid);

fid = fopen('LME_EPCA.m','w');
fprintf(fid, 'clear\n');
for i = 1:numNewLandmarks
    name = newLandmarksNames{i};

    % write s
    if size(name, 2) > 16
      fprintf(fid, '%s__s = [', name(1:16));
    else
        fprintf(fid, '%s__s = [', name);
    end
    for d = 1:dim
        fprintf(fid, '%6.32f ', s(d,i));
    end
    fprintf(fid, '];\n');

    % write M
    if size(name, 2) > 16
      fprintf(fid, '%s__M = [', name(1:16));
    else
        fprintf(fid, '%s__M = [', name);
    end
    tmp = M{i};
    for d = 1:dim
        fprintf(fid, '%6.32f ', tmp(:,d));
        fprintf(fid, ';\n');
    end
    fprintf(fid, '];\n');
end
fclose(fid);


disp('LME-EPCA coefficients have been written to disk.');
