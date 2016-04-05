% Author: Ali Ghayoor
%
% This MATLAB code is a companion for ComputeReflectiveCorrelationMetric
%
% This code reads in a csv file with 4 columns as:
% Head_angle, Bank_angle, LR_offset, reflective_correlation_metric
%
% Then, it creates a surface plot from scatter data.
% Metric value is plotted based on head an bank angles when LR offset is fixed.

M = csvread('cc_metric.csv');
LR_offset = 0;

[rows, cols] = size(M);
index([1,2,4])=true;
M0 = M(M(:,3)==LR_offset,index);
x = M0(:,1); y = M0(:,2); z = M0(:,3);

%% How do you turn a collection of XYZ triplets into a surface plot?
% The solution is to use Delaunay triangulation.
tri = delaunay(x,y);

h = trisurf(tri, x, y, z);
axis vis3d
