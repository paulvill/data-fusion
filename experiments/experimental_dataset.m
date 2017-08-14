% this script colors movies from static snapshots using semi supervised learning 
% the first step is to import and apply image preprocessing on movie frames
% and snapshots
% the scattering transformation is then computed on the entire set of
% images
% an affinity matrix is obtained from the scattering transformation and
% used to compute the semi-supervised algorithm for each channel
% finally a colored movie is stored as a sequence of images, either
% combined or channel by channel

%% Importing movie frames

% Number of movies to load.
num_movies = 7;

% Movie/snapshot loading parameters.
npixels = 512;

% Resizing for faster computation of diffusion maps.
npixels_dm = 100;

% Normalization parameters controlling the width of the averaging kernel and
% the regularization of the normalization, respectively.
sigma = 10;
epsilon = 0.01;

% The parameters for the logistic contrast boosting. Here, a is the scaling
% factor, with higher values giving a more step-like function and b is the
% offset determining which values are mapped to 0.5. Different values are
% provided for movies and snapshots.
a_movies = 10;
b_movies = 0.8;

% Specify which movies that we should use.
diffusion_movie_idx = [1 2 3 4 5 6 7]'; 

% The relative angles of the movies determined by hand.
movies_theta = [-95 -80 -80 -95 -120 -95 -55]';

% Load all the movies.
movies = empty_movie_set(npixels, 1);

fprintf('Loading movies...');
for k = 1:num_movies
    movie_opt = struct();
    movie_opt.data_dir = sprintf('data/movie%d', k);
    movie_opt.image_name = 'frame';
    movie_opt.image_ext = 'tif';
    movie_opt.dim = 2;
    movie_opt.npixels = npixels;

    movie_k = load_movie_set(movie_opt);

    movie_k.times = [1:numel(movie_k.times)]';

    movies = cat_movie_sets(movies, movie_k);
end
fprintf('OK\n');

% Save original movies for colorizing later. We don't want to
% use all the preprocessing that comes after this.
movies_orig = movies;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering movies...');
% we center movies sequence-wise to avoid jittering in the reconstruction
for i = 1:max(movies.movie_idx)
movies_orig.images(:,:,movies.movie_idx == i) = mean_center_zstack(movies_orig.images(:,:,movies.movie_idx == i),1,'TRUE');
end
fprintf('OK\n');

% Movies are rotated into reference position
fprintf('Rotating movies...');
movies_orig.images = rotate_images(movies_orig.images, ...
    movies_theta(diffusion_movie_idx(movies_orig.movie_idx)));
fprintf('OK\n');

% Resize images to make them easier to deal with computationally.
movies.images = resize_images(movies.images, npixels_dm);

% The movies and snapshots come in as 8-bit integers. For convenience, we'll
% want to treat them as floating-point numbers, however.
movies.images = double(movies.images)/256;

% To remove low-frequency variation due to illumination and defective embryos,
% we blur each image and use this to normalize the pixels of each image.
fprintf('Normalizing movies...');
movies.images = local_normalize_images(movies.images, sigma, epsilon);
fprintf('OK\n');

% To get rid of low-intensity areas and boost mid-intensity pixels, we apply
% a contrast-boosting filter in the form of a logistic function.
fprintf('Increasing contrast in movies...');
movies.images = increase_contrast_images(movies.images, a_movies, b_movies);
fprintf('OK\n');

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 1));
movies.images = matfun(center_im, movies.images, 3);
fprintf('OK\n');

% Movies are rotated into reference position
fprintf('Rotating movies...');
movies.images = rotate_images(movies.images, ...
    movies_theta(diffusion_movie_idx(movies.movie_idx)));
fprintf('OK\n');


%% Load the first set of snapshots

% Normalization parameters controlling the width of the averaging kernel and
% the regularization of the normalization, respectively.
sigma = 10;
epsilon = 0.01;

% The parameters for the logistic contrast boosting. Here, a is the scaling
% factor, with higher values giving a more step-like function and b is the
% offset determining which values are mapped to 0.5. Different values are
% provided for movies and snapshots.
a_snapshots = 10;
b_snapshots = 1.0;

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set1/';
snapshot_opt.image_name = 'ordered';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots1 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots1.times = [1:numel(snapshots1.times)]';
fprintf('OK\n');

% Save original snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig1 = snapshots1;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig1.images = matfun(center_im, snapshots_orig1.images, 4);
fprintf('OK\n');

% Resize images to make them easier to deal with computationally.
snapshots1.images = resize_images(snapshots1.images, npixels_dm);

% Extract the first channel for snapshots.
snapshots1.images = permute(snapshots1.images(:,:,1,:), [1 2 4 3]);

% The movies and snapshots come in as 8-bit integers. For convenience, we'll
% want to treat them as floating-point numbers, however.
snapshots1.images = double(snapshots1.images)/256;

% To remove low-frequency variation due to illumination and defective embryos,
% we blur each image and use this to normalize the pixels of each image.
fprintf('Normalizing images...');
snapshots1.images = local_normalize_images(snapshots1.images, sigma, epsilon);
fprintf('OK\n');

% To get rid of low-intensity areas and boost mid-intensity pixels, we apply
% a contrast-boosting filter in the form of a logistic function.
fprintf('Increasing contrast in images...');
snapshots1.images = ...
    increase_contrast_images(snapshots1.images, a_snapshots, b_snapshots);
fprintf('OK\n');

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 1));
snapshots1.images = matfun(center_im, snapshots1.images, 3);
fprintf('OK\n');



%% Load the second set of snapshots

% Normalization parameters controlling the width of the averaging kernel and
% the regularization of the normalization, respectively.
sigma = 10;
epsilon = 0.01;

% The parameters for the logistic contrast boosting. Here, a is the scaling
% factor, with higher values giving a more step-like function and b is the
% offset determining which values are mapped to 0.5. Different values are
% provided for movies and snapshots.
a_snapshots = 10;
b_snapshots = 1.0;

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set2/';
snapshot_opt.image_name = 'emb';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots2 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots2.times = [1:numel(snapshots2.times)]';
fprintf('OK\n');

snapshots2.images(:,:,[1,2,3],:) = snapshots2.images(:,:,[3,2,1],:);

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig2 = snapshots2;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig2.images = matfun(center_im, snapshots_orig2.images, 4);
fprintf('OK\n');

% Resize images to make them easier to deal with computationally.
snapshots2.images = resize_images(snapshots2.images, npixels_dm);

% Extract the first channel for both movies and snapshots.
snapshots2.images = permute(snapshots2.images(:,:,1,:), [1 2 4 3]);

% The movies and snapshots come in as 8-bit integers. For convenience, we'll
% want to treat them as floating-point numbers, however.
snapshots2.images = double(snapshots2.images)/256;

% To remove low-frequency variation due to illumination and defective embryos,
% we blur each image and use this to normalize the pixels of each image.
fprintf('Normalizing images...');
snapshots2.images = local_normalize_images(snapshots2.images, sigma, epsilon);
fprintf('OK\n');

% To get rid of low-intensity areas and boost mid-intensity pixels, we apply
% a contrast-boosting filter in the form of a logistic function.
fprintf('Increasing contrast in images...');
snapshots2.images = ...
    increase_contrast_images(snapshots2.images, a_snapshots, b_snapshots);
fprintf('OK\n');

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 1));
snapshots2.images = matfun(center_im, snapshots2.images, 3);
fprintf('OK\n');


%%% loading the other channels

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set2/';
snapshot_opt.image_name = 'ind';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots2_1 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots2_1.times = [1:numel(snapshots2_1.times)]';
fprintf('OK\n');

snapshots2_1.images(:,:,[1,2,3],:) = snapshots2_1.images(:,:,[3,2,1],:);

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig2_1 = snapshots2_1;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig2_1.images = matfun(center_im, snapshots_orig2_1.images, 4);
fprintf('OK\n');




%%  Load the third set of snapshots

% Normalization parameters controlling the width of the averaging kernel and
% the regularization of the normalization, respectively.
sigma = 5;
epsilon = 0.05;

% The parameters for the logistic contrast boosting. Here, a is the scaling
% factor, with higher values giving a more step-like function and b is the
% offset determining which values are mapped to 0.5. Different values are
% provided for movies and snapshots.
a_snapshots = 15;
b_snapshots = 0.5;


% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set3/';
snapshot_opt.image_name = 'emb';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots3 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots3.times = [1:numel(snapshots3.times)]';
fprintf('OK\n');

snapshots3.images(:,:,[1,2,3],:) = snapshots3.images(:,:,[3,2,1],:);

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig3 = snapshots3;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 0));
snapshots_orig3.images = matfun(center_im, snapshots_orig3.images, 4);
fprintf('OK\n');

% Resize images to make them easier to deal with computationally.
snapshots3.images = resize_images(snapshots3.images, npixels_dm);

% Extract the first channel for both movies and snapshots.
snapshots3.images = permute(snapshots3.images(:,:,1,:), [1 2 4 3]);

% The movies and snapshots come in as 8-bit integers. For convenience, we'll
% want to treat them as floating-point numbers, however.
snapshots3.images = double(snapshots3.images)/256;

% To remove low-frequency variation due to illumination and defective embryos,
% we blur each image and use this to normalize the pixels of each image.
fprintf('Normalizing images...');
snapshots3.images = local_normalize_images(snapshots3.images, sigma, epsilon);
fprintf('OK\n');

% To get rid of low-intensity areas and boost mid-intensity pixels, we apply
% a contrast-boosting filter in the form of a logistic function.
fprintf('Increasing contrast in images...');
snapshots3.images = ...
    increase_contrast_images(snapshots3.images, a_snapshots, b_snapshots);
fprintf('OK\n');

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 0));
snapshots3.images = matfun(center_im, snapshots3.images, 3);
fprintf('OK\n');

%%% loading the other channels

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set3/';
snapshot_opt.image_name = 'ind';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots3_1 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots3_1.times = [1:numel(snapshots3_1.times)]';
fprintf('OK\n');

snapshots3_1.images(:,:,[1,2,3],:) = snapshots3_1.images(:,:,[3,2,1],:);

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig3_1 = snapshots3_1;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig3_1.images = matfun(center_im, snapshots_orig3_1.images, 4);
fprintf('OK\n');

%% Load the fourth set of snapshots

% Normalization parameters controlling the width of the averaging kernel and
% the regularization of the normalization, respectively.
sigma = 10;
epsilon = 0.01;

% The parameters for the logistic contrast boosting. Here, a is the scaling
% factor, with higher values giving a more step-like function and b is the
% offset determining which values are mapped to 0.5. Different values are
% provided for movies and snapshots.
a_snapshots = 10;
b_snapshots = 1.0;

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set4/';
snapshot_opt.image_name = 'emb';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots4 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots4.times = [1:numel(snapshots4.times)]';
fprintf('OK\n');

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig4 = snapshots4;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig4.images = matfun(center_im, snapshots_orig4.images, 4);
fprintf('OK\n');

% Resize images to make them easier to deal with computationally.
snapshots4.images = resize_images(snapshots4.images, npixels_dm);

% Extract the first channel for both movies and snapshots.
snapshots4.images = permute(snapshots4.images(:,:,1,:), [1 2 4 3]);

% The movies and snapshots come in as 8-bit integers. For convenience, we'll
% want to treat them as floating-point numbers, however.
snapshots4.images = double(snapshots4.images)/256;

% To remove low-frequency variation due to illumination and defective embryos,
% we blur each image and use this to normalize the pixels of each image.
fprintf('Normalizing images...');
snapshots4.images = local_normalize_images(snapshots4.images, sigma, epsilon);
fprintf('OK\n');

% To get rid of low-intensity areas and boost mid-intensity pixels, we apply
% a contrast-boosting filter in the form of a logistic function.
fprintf('Increasing contrast in images...');
snapshots4.images = ...
    increase_contrast_images(snapshots4.images, a_snapshots, b_snapshots);
fprintf('OK\n');

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 1));
snapshots4.images = matfun(center_im, snapshots4.images, 3);
fprintf('OK\n');


%%% loading the other channels

% Load the snapshots.
snapshot_opt.data_dir = 'data/data_set4/';
snapshot_opt.image_name = 'ind_rho';
snapshot_opt.angles_name = 'fixed_thetas.txt';
snapshot_opt.npixels = npixels;

fprintf('Loading snapshots...');
snapshots4_1 = load_movie_set(snapshot_opt);
% Put in fake times.
snapshots4_1.times = [1:numel(snapshots4_1.times)]';
fprintf('OK\n');

% Save original movies and snapshots for colorizing later. We don't want to
% use all the preprocessing that comes after this.
snapshots_orig4_1 = snapshots4_1;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering images...');
center_im = @(im)(mean_center_image(im, 1, true, 15));
snapshots_orig4_1.images = matfun(center_im, snapshots_orig4_1.images, 4);
fprintf('OK\n');

%% In this part of the code, we compute the scatter transformation and the 
% affinity matrix that will then be used for semi-supervised learning

% We first remove some frames of the movies that are outlying and restricting the
% window of observation from shortly before gastrulation to shortly after 
% gastrulation to avoid complicated morphological changes which are not
% well captured in cross section and too static morphology before
% gastrulation

mask = [];

indtemp1 = find(movies.movie_idx == 1);
indtemp1 = indtemp1(90:149);
indtemp1(41:46) = [];
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 2);
indtemp1 = indtemp1(1:49);
indtemp1(31:41)= [];
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 3);
indtemp1 = indtemp1(13:44);
indtemp1(12:26) = [];
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 4);
indtemp1 = indtemp1(1:43);
indtemp1(30:34) = [];
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 5);
indtemp1 = indtemp1(1:42);
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 6);
indtemp1 = indtemp1(1:58);
mask = [mask;indtemp1];

indtemp1 = find(movies.movie_idx == 7);
indtemp1 = indtemp1(1:63);
indtemp1(17) = [];
mask = [mask;indtemp1];

% the outliers are removed
movies_sub = subset_movie_set(movies, mask);
movies_orig_sub = subset_movie_set(movies_orig, mask);

% Put the movies and snapshots together into one movie set to simplify later
% calculations.
all = cat_movie_sets(movies_sub, snapshots1);
all = cat_movie_sets(all, snapshots2);
all = cat_movie_sets(all, snapshots3);
all = cat_movie_sets(all, snapshots4);

% The order (number of layers) and scale (2^J is the averaging window size)
% of the scattering transform. For M = 0, only blurred images are given,
% while M = 1 and M = 2 retains more information on the finer-scale spatial
% structure of the images.
scat_M = 1;
scat_J = 6;

% To obtain some invariance to translation (our centering above may not be as
% good as we need it to be) and create stability to deformation (in particular
% some snapshots can be very deformed from the ideal circular embryo), we
% compute the scattering transform of the images.
fprintf('Calculating scattering transforms...');
pad_factor = 2;
orig_sz = size(all.images(:,:,1));
padded_sz = pad_factor*orig_sz;
scat_opt.M = scat_M;
filt_opt.J = scat_J;
Wop = wavelet_factory_2d(padded_sz, filt_opt, scat_opt);

S_all = scat_images(all.images, Wop, false, pad_factor);
S_all = reshape(S_all, [], size(S_all, 4));
fprintf('OK\n');


% Center each of datasets to have mean zero.
fprintf('Centering feature vectors...');
for i = 1:length(unique(all.movie_idx)),
    indtemp = find(all.movie_idx == i);
    S_all(:,indtemp) = bsxfun(@minus,S_all(:,indtemp), mean(S_all(:,indtemp),2));
end
fprintf('OK\n');

% Calculate the pairwise distances
fprintf('Calculating pairwise differences ...');
V = calc_pairwise_distances(S_all);
fprintf('OK\n');

min_distance = min(V+diag(Inf(size(V, 1), 1)));
median_min_distance = median(min_distance);

mask = min_distance<3*median_min_distance;
V = V(mask,mask);

clear embed_coords;
fprintf('Calculating diffusion maps for sanity check ...');
[~, embed_coords] = DiffusionMapsFromDistanceGlobal(V, 1, 50);
embed_coords = embed_coords(:,2:4);
fprintf('OK\n');

% figure showing low dimensional embedding of the cloud of points in the
% diffusion maps space
figure, 
subplot(2,2,1)
scatter3(embed_coords(:,1),embed_coords(:,2),embed_coords(:,3),50,all.movie_idx(mask),'filled');
grid off
xlabel('Diff. Coord. 1')
ylabel('Diff. Coord. 2')
zlabel('Diff. Coord. 3')
subplot(2,2,2)
scatter(embed_coords(:,1),embed_coords(:,2),50,all.movie_idx(mask),'filled');
grid off
xlabel('Diff. Coord. 1')
ylabel('Diff. Coord. 2')
subplot(2,2,3)
scatter(embed_coords(:,1),embed_coords(:,3),50,all.movie_idx(mask),'filled');
grid off
xlabel('Diff. Coord. 1')
ylabel('Diff. Coord. 3')
subplot(2,2,4)
scatter(embed_coords(:,2),embed_coords(:,3),50,all.movie_idx(mask),'filled');
grid off
xlabel('Diff. Coord. 2')
ylabel('Diff. Coord. 3')


%% In this section, we visualize the pairwise differences and the affinity matrix

% Calculate the pairwise distances 
fprintf('Calculating pairwise differences ...');
V_vis = calc_pairwise_distances(S_all);
fprintf('OK\n');

% figure showing the pairwise differences between images
figure, 
imagesc(V_vis);
hold on,
for i = 1:max(all.movie_idx),
    indtemp = find(all.movie_idx == i);    
    rectangle('Position',[indtemp(1) indtemp(1) indtemp(end)-indtemp(1) indtemp(end)-indtemp(1)],'EdgeColor','w','LineWidth',4)
end

% figure showing the affinity matrix
W_vis = AffinityFromDistance(V_vis,10);
figure, 
imagesc(W_vis);
hold on,
for i = 1:max(all.movie_idx),
    indtemp = find(all.movie_idx == i);    
    rectangle('Position',[indtemp(1) indtemp(1) indtemp(end)-indtemp(1) indtemp(end)-indtemp(1)],'EdgeColor','w','LineWidth',4)
end

%% Predicting dpERK on the unlabeled movie frames

t = cputime;
ntot = sum(mask);
nmov = length(find(find(all.movie_idx(mask) < 8)));
W = AffinityFromDistance(V,10);

indl = find(all.movie_idx(mask) == 8 | all.movie_idx(mask) == 9 | all.movie_idx(mask) == 10)';
indu = find(all.movie_idx(mask) ~= 8 & all.movie_idx(mask) ~= 9 & all.movie_idx(mask) ~= 10)';
l1 = length(indl);
l2 = length(indu);

[inv_u] = ssl_estimate( W, indl, indu);

npixels = 512;
movies_colored_dpERK = movies_orig_sub;
movies_colored_dpERK.images = permute(movies_colored_dpERK.images, [1 2 4 3]);
movies_colored_dpERK.images = cat(3, movies_colored_dpERK.images, zeros(size(movies_colored_dpERK.images)));

e = cputime - t;

fprintf(['Initialization OK in ',num2str(e),' s \n']);

snapshots1_mask = mask(all.movie_idx == 8);
snapshots2_mask = mask(all.movie_idx == 9);
snapshots3_mask = mask(all.movie_idx == 10);
snapshots4_mask = mask(all.movie_idx == 11);

Y1 = cat(4, snapshots_orig1.images(:,:,2,snapshots1_mask), ...
    snapshots_orig2_1.images(:,:,3,snapshots2_mask), ...
    snapshots_orig3_1.images(:,:,3,snapshots3_mask));
Y1 = double(Y1);
Y1 = permute(Y1, [4 1 2 3]);
[Y1, sz_roll] = unroll_dim(Y1, 2);

fu1 = inv_u*Y1;
fu1 = fu1(1:nmov,:);
fu1 = roll_dim(fu1, sz_roll);
fu1 = permute(fu1, [2 3 1]);

movies_colored_dpERK.images(:,:,2,all.movie_idx<8&mask') = fu1;

e = cputime - t;

fprintf(['Movies colored in ',num2str(e),' s \n']);

movies_colored_dpERK.images = double(movies_colored_dpERK.images);

for movie_id = unique(movies_colored_dpERK.movie_idx)'
    movie_idx = find(movies_colored_dpERK.movie_idx==movie_id);

    mean_intensities = mean(mean(mean(movies_colored_dpERK.images(:,:,:,movie_idx), 1), 2), 4);

    movies_colored_dpERK.images(:,:,1,movie_idx) = 1./mean_intensities(1)*movies_colored_dpERK.images(:,:,1,movie_idx);
    movies_colored_dpERK.images(:,:,2,movie_idx) = 1./mean_intensities(2)*movies_colored_dpERK.images(:,:,2,movie_idx);
end

movies_colored_dpERK.images = movies_colored_dpERK.images/max(movies_colored_dpERK.images(:));


%% Predicting Twi on the unlabeled movie frames

t = cputime;
ntot = sum(mask);
nmov = length(find(find(all.movie_idx(mask) < 8)));
W = AffinityFromDistance(V,10);

indl = find(all.movie_idx(mask) == 8 | all.movie_idx(mask) == 11)';
indu = find(all.movie_idx(mask) ~= 8 & all.movie_idx(mask) ~= 11)';

l1 = length(indl);
l2 = length(indu);

[inv_u] = ssl_estimate( W, indl, indu);

npixels = 512;
movies_colored_twi = movies_orig_sub;
movies_colored_twi.images = permute(movies_colored_twi.images, [1 2 4 3]);
movies_colored_twi.images = cat(3, movies_colored_twi.images, zeros(size(movies_colored_twi.images)));
im_colored = zeros(npixels, npixels, 1);

e = cputime - t;

fprintf(['Initialization OK in ',num2str(e),' s \n']);

snapshots1_mask = mask(all.movie_idx == 8);
snapshots2_mask = mask(all.movie_idx == 9);
snapshots3_mask = mask(all.movie_idx == 10);
snapshots4_mask = mask(all.movie_idx == 11);

Y1 = cat(4, snapshots_orig1.images(:,:,3,snapshots1_mask), ...
    snapshots_orig4_1.images(:,:,3,snapshots4_mask));
Y1 = double(Y1);
Y1 = permute(Y1, [4 1 2 3]);
[Y1, sz_roll] = unroll_dim(Y1, 2);

fu1 = inv_u*Y1;
fu1 = fu1(1:nmov,:);
fu1 = roll_dim(fu1, sz_roll);
fu1 = permute(fu1, [2 3 1]);

movies_colored_twi.images(:,:,2,all.movie_idx<8&mask') = fu1;

e = cputime - t;

fprintf(['Movies colored in ',num2str(e),' s \n']);

movies_colored_twi.images = double(movies_colored_twi.images);

for movie_id = unique(movies_colored_twi.movie_idx)'
    movie_idx = find(movies_colored_twi.movie_idx==movie_id);

    mean_intensities = mean(mean(mean(movies_colored_twi.images(:,:,:,movie_idx), 1), 2), 4);

    movies_colored_twi.images(:,:,1,movie_idx) = 1./mean_intensities(1)*movies_colored_twi.images(:,:,1,movie_idx);
    movies_colored_twi.images(:,:,2,movie_idx) = 1./mean_intensities(2)*movies_colored_twi.images(:,:,2,movie_idx);
end

movies_colored_twi.images = movies_colored_twi.images/max(movies_colored_twi.images(:));


%% Predicting ind on the unlabeled movie frames

t = cputime;
ntot = sum(mask);
nmov = length(find(find(all.movie_idx(mask) < 8)));
W = AffinityFromDistance(V,10);

indl = find(all.movie_idx(mask) == 9 | all.movie_idx(mask) == 11)';
indu = find(all.movie_idx(mask) ~= 9 & all.movie_idx(mask) ~= 11)';

l1 = length(indl);
l2 = length(indu);

[inv_u] = ssl_estimate( W, indl, indu);

npixels = 512;
movies_colored_ind = movies_orig_sub;
movies_colored_ind.images = permute(movies_colored_ind.images, [1 2 4 3]);
movies_colored_ind.images = cat(3, movies_colored_ind.images, zeros(size(movies_colored_ind.images)));
im_colored = zeros(npixels, npixels, 1);

e = cputime - t;

fprintf(['Initialization OK in ',num2str(e),' s \n']);

snapshots1_mask = mask(all.movie_idx == 8);
snapshots2_mask = mask(all.movie_idx == 9);
snapshots3_mask = mask(all.movie_idx == 10);
snapshots4_mask = mask(all.movie_idx == 11);

Y1 = cat(4, snapshots_orig2.images(:,:,2,snapshots2_mask), ...
    snapshots_orig4_1.images(:,:,1,snapshots4_mask));
Y1 = double(Y1);
Y1 = permute(Y1, [4 1 2 3]);
[Y1, sz_roll] = unroll_dim(Y1, 2);

fu1 = inv_u*Y1;
fu1 = fu1(1:nmov,:);
fu1 = roll_dim(fu1, sz_roll);
fu1 = permute(fu1, [2 3 1]);

movies_colored_ind.images(:,:,2,all.movie_idx<8&mask') = fu1;

e = cputime - t;

fprintf(['Movies colored in ',num2str(e),' s \n']);

movies_colored_ind.images = double(movies_colored_ind.images);

for movie_id = unique(movies_colored_ind.movie_idx)'
    movie_idx = find(movies_colored_ind.movie_idx==movie_id);

    mean_intensities = mean(mean(mean(movies_colored_ind.images(:,:,:,movie_idx), 1), 2), 4);

    movies_colored_ind.images(:,:,1,movie_idx) = 1./mean_intensities(1)*movies_colored_ind.images(:,:,1,movie_idx);
    movies_colored_ind.images(:,:,2,movie_idx) = 1./mean_intensities(2)*movies_colored_ind.images(:,:,2,movie_idx);
end

movies_colored_ind.images = movies_colored_ind.images/max(movies_colored_ind.images(:));


%% Predicting Dorsal (Dl) on the unlabeled movie frames

t = cputime;
ntot = sum(mask);
nmov = length(find(find(all.movie_idx(mask) < 8)));
W = AffinityFromDistance(V,10);

indl = find(all.movie_idx(mask) == 9)';
indu = find(all.movie_idx(mask) ~= 9)';

l1 = length(indl);
l2 = length(indu);

[inv_u] = ssl_estimate( W, indl, indu);

npixels = 512;
movies_colored_dl = movies_orig_sub;
movies_colored_dl.images = permute(movies_colored_dl.images, [1 2 4 3]);
movies_colored_dl.images = cat(3, movies_colored_dl.images, zeros(size(movies_colored_dl.images)));
im_colored = zeros(npixels, npixels, 1);

e = cputime - t;

fprintf(['Initialization OK in ',num2str(e),' s \n']);

snapshots1_mask = mask(all.movie_idx == 8);
snapshots2_mask = mask(all.movie_idx == 9);
snapshots3_mask = mask(all.movie_idx == 10);
snapshots4_mask = mask(all.movie_idx == 11);

Y1 = snapshots_orig2_1.images(:,:,1,snapshots2_mask);
Y1 = double(Y1);
Y1 = permute(Y1, [4 1 2 3]);
[Y1, sz_roll] = unroll_dim(Y1, 2);

fu1 = inv_u*Y1;
fu1 = fu1(1:nmov,:);
fu1 = roll_dim(fu1, sz_roll);
fu1 = permute(fu1, [2 3 1]);

movies_colored_dl.images(:,:,2,all.movie_idx<8&mask') = fu1;

e = cputime - t;

fprintf(['Movies colored in ',num2str(e),' s \n']);

movies_colored_dl.images = double(movies_colored_dl.images);

for movie_id = unique(movies_colored_dl.movie_idx)'
    movie_idx = find(movies_colored_dl.movie_idx==movie_id);

    mean_intensities = mean(mean(mean(movies_colored_dl.images(:,:,:,movie_idx), 1), 2), 4);

    movies_colored_dl.images(:,:,1,movie_idx) = 1./mean_intensities(1)*movies_colored_dl.images(:,:,1,movie_idx);
    movies_colored_dl.images(:,:,2,movie_idx) = 1./mean_intensities(2)*movies_colored_dl.images(:,:,2,movie_idx);
end

movies_colored_dl.images = movies_colored_dl.images/max(movies_colored_dl.images(:));

%% Predicting rho on the unlabeled movie frames

t = cputime;
ntot = sum(mask);
nmov = length(find(find(all.movie_idx(mask) < 8)));
W = AffinityFromDistance(V,10);

indl = find(all.movie_idx(mask) >9)';
indu = find(all.movie_idx(mask) < 10)';

l1 = length(indl);
l2 = length(indu);


[inv_u] = ssl_estimate( W, indl, indu);

npixels = 512;
movies_colored_rho = movies_orig_sub;
movies_colored_rho.images = permute(movies_colored_rho.images, [1 2 4 3]);
movies_colored_rho.images = cat(3, movies_colored_rho.images, zeros(size(movies_colored_rho.images)));
im_colored = zeros(npixels, npixels, 1);

e = cputime - t;

fprintf(['Initialization OK in ',num2str(e),' s \n']);

snapshots1_mask = mask(all.movie_idx == 8);
snapshots2_mask = mask(all.movie_idx == 9);
snapshots3_mask = mask(all.movie_idx == 10);
snapshots4_mask = mask(all.movie_idx == 11);

Y1 = cat(4, snapshots_orig3_1.images(:,:,1,snapshots3_mask), snapshots_orig4_1.images(:,:,2,snapshots4_mask));
Y1 = double(Y1);
Y1 = permute(Y1, [4 1 2 3]);
[Y1, sz_roll] = unroll_dim(Y1, 2);

fu1 = inv_u*Y1;
fu1 = fu1(1:nmov,:);
fu1 = roll_dim(fu1, sz_roll);
fu1 = permute(fu1, [2 3 1]);

movies_colored_rho.images(:,:,2,all.movie_idx<8&mask') = fu1;

e = cputime - t;

fprintf(['Movies colored in ',num2str(e),' s \n']);

movies_colored_rho.images = double(movies_colored_rho.images);

for movie_id = unique(movies_colored_rho.movie_idx)'
    movie_idx = find(movies_colored_rho.movie_idx==movie_id);

    mean_intensities = mean(mean(mean(movies_colored_rho.images(:,:,:,movie_idx), 1), 2), 4);

    movies_colored_rho.images(:,:,1,movie_idx) = 1./mean_intensities(1)*movies_colored_rho.images(:,:,1,movie_idx);
    movies_colored_rho.images(:,:,2,movie_idx) = 1./mean_intensities(2)*movies_colored_rho.images(:,:,2,movie_idx);
end

movies_colored_rho.images = movies_colored_rho.images/max(movies_colored_rho.images(:));

%% Coloring the movie frames according to the color code

% color code
% dpERK = [213, 36, 2];
% Twi = [92, 188, 146];
% Ind = [103, 139, 191];
% Dorsal = [229, 195, 207];
% Rho = [255, 223, 75];
col = [[213, 36, 2];[92, 188, 146];[103, 139, 191];[229, 195, 207];[255, 223, 75]]/255;

nuclei = zeros(size(movies_colored_dpERK.images(:,:,:,:)));
nuclei(:,:,1,:) = movies_colored_dpERK.images(:,:,1,:);
nuclei(:,:,2,:) = movies_colored_dpERK.images(:,:,1,:);
nuclei(:,:,3,:) = movies_colored_dpERK.images(:,:,1,:);

dpERK = zeros(size(movies_colored_dpERK.images(:,:,:,:)));
dpERK(:,:,1,:) = col(1,1)*movies_colored_dpERK.images(:,:,2,:);
dpERK(:,:,2,:) = col(1,2)*movies_colored_dpERK.images(:,:,2,:);
dpERK(:,:,3,:) = col(1,3)*movies_colored_dpERK.images(:,:,2,:);

twi = zeros(size(movies_colored_twi.images(:,:,:,:)));
twi(:,:,1,:) = col(2,1)*movies_colored_twi.images(:,:,2,:);
twi(:,:,2,:) = col(2,2)*movies_colored_twi.images(:,:,2,:);
twi(:,:,3,:) = col(2,3)*movies_colored_twi.images(:,:,2,:);

ind = zeros(size(movies_colored_ind.images(:,:,:,:)));
ind(:,:,1,:) = col(3,1)*movies_colored_ind.images(:,:,2,:);
ind(:,:,2,:) = col(3,2)*movies_colored_ind.images(:,:,2,:);
ind(:,:,3,:) = col(3,3)*movies_colored_ind.images(:,:,2,:);

dl = zeros(size(movies_colored_dl.images(:,:,:,:)));
dl(:,:,1,:) = col(4,1)*movies_colored_dl.images(:,:,2,:);
dl(:,:,2,:) = col(4,2)*movies_colored_dl.images(:,:,2,:);
dl(:,:,3,:) = col(4,3)*movies_colored_dl.images(:,:,2,:);

rho = zeros(size(movies_colored_rho.images(:,:,:,:)));
rho(:,:,1,:) = col(5,1)*movies_colored_rho.images(:,:,2,:);
rho(:,:,2,:) = col(5,2)*movies_colored_rho.images(:,:,2,:);
rho(:,:,3,:) = col(5,3)*movies_colored_rho.images(:,:,2,:);


fused = [];

fused = nuclei(:,:,:,:) + twi(:,:,:,:) + ind(:,:,:,:) + dl(:,:,:,:) + rho(:,:,:,:) + dpERK(:,:,:,:);
fused(:,:,1,:) = fused(:,:,1,:)/max(max(max(fused(:,:,1,:))));
fused(:,:,2,:) = fused(:,:,2,:)/max(max(max(fused(:,:,2,:))));
fused(:,:,3,:) = fused(:,:,3,:)/max(max(max(fused(:,:,3,:))));

%% save the fused movie

%  with all the channels together
d = fullfile('output', 'coloredmovie');
mkdirp(d);

indtemp = find(all.movie_idx == 5);
figure, 
for k = 1:length(indtemp),
    imagesc(fused(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(fused(:,:,:,indtemp(k)),fullfile(d, sprintf('%04d.png', k)));
end

% save the fused movie with all channels separated as gray scale

d = fullfile('output', 'grayscale');
mkdirp(d);

indtemp = find(all.movie_idx == 5);

% nuclei
imtemp = [];
imtemp = movies_colored_dpERK.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_dpERK.images(:,:,1,:);
imtemp(:,:,2,:) = movies_colored_dpERK.images(:,:,1,:);
imtemp(:,:,3,:) = movies_colored_dpERK.images(:,:,1,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('nuclei%04d.png', k)));
end

% dpERK

imtemp = [];
imtemp = movies_colored_dpERK.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_dpERK.images(:,:,2,:);
imtemp(:,:,2,:) = movies_colored_dpERK.images(:,:,2,:);
imtemp(:,:,3,:) = movies_colored_dpERK.images(:,:,2,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('dpERK%04d.png', k)));
end

% twist 

imtemp = [];
imtemp = movies_colored_twi.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_twi.images(:,:,2,:);
imtemp(:,:,2,:) = movies_colored_twi.images(:,:,2,:);
imtemp(:,:,3,:) = movies_colored_twi.images(:,:,2,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('twist%04d.png', k)));
end

% ind

imtemp = [];
imtemp = movies_colored_ind.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_ind.images(:,:,2,:);
imtemp(:,:,2,:) = movies_colored_ind.images(:,:,2,:);
imtemp(:,:,3,:) = movies_colored_ind.images(:,:,2,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('ind%04d.png', k)));
end

% dorsal

imtemp = [];
imtemp = movies_colored_dl.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_dl.images(:,:,2,:);
imtemp(:,:,2,:) = movies_colored_dl.images(:,:,2,:);
imtemp(:,:,3,:) = movies_colored_dl.images(:,:,2,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('dorsal%04d.png', k)));
end

% rhomboid

imtemp = [];
imtemp = movies_colored_rho.images(:,:,:,:);
imtemp(:,:,1,:) = movies_colored_rho.images(:,:,2,:);
imtemp(:,:,2,:) = movies_colored_rho.images(:,:,2,:);
imtemp(:,:,3,:) = movies_colored_rho.images(:,:,2,:);

figure, 
for k = 1:length(indtemp),
    imagesc(imtemp(:,:,:,indtemp(k)));
    pause(0.1);
    imwrite(imtemp(:,:,:,indtemp(k)),fullfile(d, sprintf('rhomboid%04d.png', k)));
end


