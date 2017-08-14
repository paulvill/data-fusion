% Run toy examples in the `examples` folder.

% Apply our data fusion method to the synthetic spiral dataset, where a one-
% dimensional trajectory embedded in three dimensions is reconstructed from
% its projection onto two dimensions.
spiral_example;

% Perform the same experiment, but for different numbers numbers of unlabeled
% data points. The errors are cross-validated over the entire dataset, providing
% a measure of variability in the reconstruction.
spiral_example_cross_validation;
