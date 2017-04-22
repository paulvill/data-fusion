function PCA_data = create_PCA_data(images)

[npixels, npixels2, nimages] = size(images);

PCA_data = reshape(double(images), npixels*npixels2, nimages)';