root_dir = fileparts(mfilename('fullpath'));

% Subfolders
addpath(fullfile(root_dir, 'image'));
addpath(fullfile(root_dir, 'utilities'));
addpath(fullfile(root_dir, 'examples'));
addpath(fullfile(root_dir, 'experiments'));
addpath(fullfile(root_dir, 'manifold'));
addpath(fullfile(root_dir, 'visualization'));

% Scripts from Lim, B. et al., Current Biology, 2015
addpath(fullfile(root_dir, 'DLL15'));
% Scripts from Lederman, R.R. and Talmon, R., Tech. Rep., Yale, 2014
addpath(fullfile(root_dir, 'LT14'));

% ScatNet Toolbox for Scattering Transform
addpath(fullfile(root_dir, 'scatnet-0.2'));
addpath_scatnet;

clear root_dir;
