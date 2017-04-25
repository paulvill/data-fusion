root_dir = fileparts(mfilename('fullpath'));

% Scripts from Lim, B. et al., Current Biology, 2015
addpath(fullfile(root_dir, 'DLL15'));
% Scripts from Dsilva, C.J. et al., Development, 2015
addpath(fullfile(root_dir, 'LDL15'));
% Scripts from Lederman, R.R. and Talmon, R., Tech. Rep., Yale, 2014
addpath(fullfile(root_dir, 'LT14'));

% ScatNet Toolbox for Scattering Transform
addpath(fullfile(root_dir, 'scatnet-0.2'));
addpath_scatnet;

% Visualization settings
presets;
