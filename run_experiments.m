% Run scripts for experimental data in the `experiments` folder.

% First, download the dataset.
download_data;

% Apply our data fusion method to the dataset of experimental images depicting
% fly embryos at different stages of gastrulation. Here, the mapping between
% nuclear morphology to spatial distribution of different chemical species is
% learned from snapshots and used to color a live imaging movie.
experimental_dataset;

% The same experiment is performed, but for different subsets of the dataset,
% characterizing the variability in the reconstruction due to the sampling of
% the datset.
experimental_dataset_cross_validation;
