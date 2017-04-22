addpath_datafusion

% this script aims at computing the K-fold cross-validation error for each
% of the labels, for each of the datasets
% we first import all the images, movies and snapshots and apply image
% preprocessing steps on them
% we then compute the scatter transformation of the images and the affinity
% matrix
% we finally compute the cross validation normalized absolute error for
% each of the channels and each of the datasets, with a varying number of
% unlabeled movie frames

%% Importing movie frames

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
movie_opt.data_dir = 'data/movies';
movie_opt.image_ext = 'avi';
movie_opt.dim = 3;
movie_opt.npixels = npixels;

fprintf('Loading movies...');
movies = load_movie_set(movie_opt);
fprintf('OK\n');

% Suppress 2nd and 3rd channels. These have no real information.
movies.images(:,:,2:3,:) = 0;

% Save original movies for colorizing later. We don't want to
% use all the preprocessing that comes after this.
movies_orig = movies;

% To make sure all images are aligned prior to comparing them, we center them.
fprintf('Centering movies...');
% we center movies sequence-wise to avoid jittering in the reconstruction
for i = 1:max(movies.movie_idx)
movies_orig.images(:,:,:,movies.movie_idx == i) = mean_center_zstack(movies_orig.images(:,:,:,movies.movie_idx == i),1,'TRUE');
end
fprintf('OK\n');

% Movies are rotated into reference position
fprintf('Rotating movies...');
movies_orig.images = rotate_images(movies_orig.images, ...
    movies_theta(diffusion_movie_idx(movies_orig.movie_idx)));
fprintf('OK\n');

movies.images_orig = movies_orig.images;

% Resize images to make them easier to deal with computationally.
movies.images = resize_images(movies.images, npixels_dm);

% Extract the first channel for movies
movies.images = permute(movies.images(:,:,1,:), [1 2 4 3]);

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

snapshots1.images_orig = snapshots_orig1.images;

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

snapshots2.images_orig = snapshots_orig2.images;

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

snapshots3.images(:,:,[1,2,3],:) = snapshots3.images(:,:,[3,2,1],:)

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

snapshots3_1.images(:,:,[1,2,3],:) = snapshots3_1.images(:,:,[3,2,1],:)

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

snapshots4.images_orig = snapshots_orig4.images;

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

%% Here we set up the parameters of cross-validation, depending on the 
% number of unlabeled data points and the number of repetition, this part
% can take some time

% number of unlabeled data points
n_unlbds = [(0:10:50),(100:50:250),length(find(all.movie_idx<8))];
M = length(n_unlbds);

% number of repetitions
nrep = 10;

% the normalized absolute error is computed for each 2 or 3 channels of each 4
% datasets
abs_error = zeros(M,nrep,3,4);


%% K-fold cross validation and varying the number of unlabeled samples - data set 1

K = 6;

% varying the number of unlabeled samples
for u = 1:M
    fprintf([num2str(n_unlbds(u)),' \n']);
    
    % repetitions over the process of recoloring with K bins randomly drawn
    % and a fixed number of unlabeled samples
    for r = 1:nrep
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = find(all.movie_idx ~= 8 )';
        indl = find(all.movie_idx == 8 )';
        
        mask_unlbld = randperm(length(indu),nb_unlbld);
        indu_mask = sort(indu(mask_unlbld));
        mask_tot = [indu_mask,indl];
        nmov = length(indu_mask);
        
        W = AffinityFromDistance(V,10);
        
        % defining the K subsamples of the labeled data points
        nb_labels = length(snapshots_orig1.movie_idx);
        ind_l_tot = 1:nb_labels;
        
        ind_l_sub = zeros(nb_labels/K,K);
        mask_label = ind_l_tot;
        
        for k = 1:K
            ind_temp = randperm(nb_labels - (k-1)*nb_labels/K,nb_labels/K);
            ind_l_sub(:,k) = mask_label(ind_temp);
            mask_label(ind_temp) = [];
        end
        
        t = cputime;
        
        % color the selected snapshots as if they were unlabeled
        % dpERK
        snapshots_recolored1 = zeros(size(snapshots_orig1.images(:,:,2,:)));
        im_orig1 = double(snapshots_orig1.images(:,:,2,:));
        
        % Twist
        snapshots_recolored2 = zeros(size(snapshots_orig1.images(:,:,3,:)));
        im_orig2 = double(snapshots_orig1.images(:,:,3,:));
        
        % for each of the bins
        for k = 1:K
            indl_k = indl;
            
            indu_k = [indu,indl(ind_l_sub(:,k))];
            indl_k(ind_l_sub(:,k)) = [];
            
            l1 = length(indl_k);
            l2 = length(indu_k);
            
            % the semi-supervised learning problem is solved on each
            % combination of labeled and unlabeled data points
            [inv_u] = ssl_estimate( W, indl_k, indu_k);
            
            npixels = 512;
            
            e = cputime - t;
            fprintf([num2str(k),'/',num2str(K),' - Initialization - OK in ',num2str(e),' s \n']);
            t = cputime;
            
            for i = 1:npixels
                for j = 1:npixels
                    
                    Y1 = [double(squeeze(snapshots_orig1.images(i,j,2,indl_k-min(indl_k)+1)))];
                    fu1 = inv_u*Y1;
                    
                    snapshots_recolored1(i,j,ind_l_sub(:,k)) = fu1(end-length(ind_l_sub(:,k))+1:end);
                    
                    Y2 = [double(squeeze(snapshots_orig1.images(i,j,3,indl_k-min(indl_k)+1)))];
                    fu2 = inv_u*Y2;
                    
                    snapshots_recolored2(i,j,ind_l_sub(:,k)) = fu2(end-length(ind_l_sub(:,k))+1:end);
                end
            end
            
            e = cputime - t;
            fprintf(['Snapshots recolored in ',num2str(e),' s \n']);
            t = cputime;
            
        end
        
        % the recolored snapshots are compared to the original snapshots
        max_int = max(im_orig1(:));
        min_int = min(im_orig1(:));
        
        abs_error(u,r,1,1) = (1/length(im_orig1(:)))*sum(abs(snapshots_recolored1(:) - im_orig1(:)))/(max_int - min_int);
        
        max_int = max(im_orig2(:));
        min_int = min(im_orig2(:));
        
        abs_error(u,r,2,1) = (1/length(im_orig2(:)))*sum(abs(snapshots_recolored2(:) - im_orig2(:)))/(max_int - min_int);
        
    end
end

%% K-fold cross validation and varying the number of unlabeled samples - data set 2
K = 3;

% varying the number of unlabeled samples
for u = 1:M
    fprintf([num2str(n_unlbds(u)),' \n']);
    
    % repetitions over the process of recoloring with K bins randomly drawn
    % and a fixed number of unlabeled samples
    for r = 1:nrep
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = find(all.movie_idx ~= 9 )';
        indl = find(all.movie_idx == 9 )';
        
        mask_unlbld = randperm(length(indu),nb_unlbld);
        indu_mask = sort(indu(mask_unlbld));
        mask_tot = [indu_mask,indl];
        nmov = length(indu_mask);
        
        W = AffinityFromDistance(V,10);
        
        % defining the K subsamples of the labeled data points
        nb_labels = length(snapshots_orig2.movie_idx);
        ind_l_tot = 1:nb_labels;
        
        ind_l_sub = {};%zeros(nb_labels/K,K);
        mask_label = ind_l_tot;
        
        for k = 1:K-1,
            ind_temp = randperm(nb_labels - (k-1)*round(nb_labels/K),round(nb_labels/K));
            ind_l_sub{k} = mask_label(ind_temp);
            mask_label(ind_temp) = [];
        end
        ind_l_sub{K} = mask_label;
        
        t = cputime;
        % color the selected snapshots as if they were unlabeled
        % dpERK
        snapshots_recolored1 = zeros(size(snapshots_orig2_1.images(:,:,3,:)));
        im_orig1 = double(snapshots_orig2_1.images(:,:,3,:));
        
        % ind
        snapshots_recolored2 = zeros(size(snapshots_orig2_1.images(:,:,2,:)));
        im_orig2 = double(snapshots_orig2_1.images(:,:,2,:));
        
        % Dorsal
        snapshots_recolored3 = zeros(size(snapshots_orig2_1.images(:,:,1,:)));
        im_orig3 = double(snapshots_orig2_1.images(:,:,1,:));
        
        for k = 1:K
            indl_k = indl;
            
            indu_k = [indu,indl(ind_l_sub{k})];
            indl_k(ind_l_sub{k}) = [];
            
            l1 = length(indl_k);
            l2 = length(indu_k);
            
            % the semi-supervised learning problem is solved on each
            % combination of labeled and unlabeled data points
            [inv_u] = ssl_estimate( W, indl_k, indu_k);
            
            npixels = 512;
            
            e = cputime - t;
            fprintf([num2str(k),'/',num2str(K),' - Initialization - OK in ',num2str(e),' s \n']);
            t = cputime;
            
            for i = 1:npixels
                for j = 1:npixels
                    
                    Y1 = [double(squeeze(im_orig1(i,j,indl_k-min(indl_k)+1)))];
                    fu1 = inv_u*Y1;
                    
                    snapshots_recolored1(i,j,ind_l_sub{k}) = fu1(end-length(ind_l_sub{k})+1:end);
                    
                    Y2 = [double(squeeze(im_orig2(i,j,indl_k-min(indl_k)+1)))];
                    fu2 = inv_u*Y2;
                    
                    snapshots_recolored2(i,j,ind_l_sub{k}) = fu2(end-length(ind_l_sub{k})+1:end);
                    
                    Y3 = [double(squeeze(im_orig3(i,j,indl_k-min(indl_k)+1)))];
                    fu3 = inv_u*Y3;
                    
                    snapshots_recolored3(i,j,ind_l_sub{k}) = fu3(end-length(ind_l_sub{k})+1:end);
                end
            end
            
            e = cputime - t;
            fprintf(['Snapshots recolored in ',num2str(e),' s \n']);
            t = cputime;
            
        end
        
        % the recolored snapshots are compared to the original snapshots
        max_int = max(im_orig1(:));
        min_int = min(im_orig1(:));
        
        abs_error(u,r,1,2) = (1/length(im_orig1(:)))*sum(abs(snapshots_recolored1(:) - im_orig1(:)))/(max_int - min_int);
        
        max_int = max(im_orig2(:));
        min_int = min(im_orig2(:));
        
        abs_error(u,r,2,2) = (1/length(im_orig2(:)))*sum(abs(snapshots_recolored2(:) - im_orig2(:)))/(max_int - min_int);
        
        max_int = max(im_orig3(:));
        min_int = min(im_orig3(:));
        
        abs_error(u,r,3,2) = (1/length(im_orig3(:)))*sum(abs(snapshots_recolored3(:) - im_orig3(:)))/(max_int - min_int);
    end
end



%% K-fold cross validation and varying the number of unlabeled samples - data set 3

K = 3;

% varying the number of unlabeled samples
for u = [1,M]
    fprintf([num2str(n_unlbds(u)),' \n']);
    
    % repetitions over the process of recoloring with K bins randomly drawn
    % and a fixed number of unlabeled samples
    for r = 1:nrep
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = find(all.movie_idx ~= 10 )';
        indl = find(all.movie_idx == 10 )';
        
        mask_unlbld = randperm(length(indu),nb_unlbld);
        indu_mask = sort(indu(mask_unlbld));
        mask_tot = [indu_mask,indl];
        nmov = length(indu_mask);
        
        W = AffinityFromDistance(V,10);
        
        % defining the K subsamples of the labeled data points
        nb_labels = length(snapshots_orig3.movie_idx);
        ind_l_tot = 1:nb_labels;
        
        ind_l_sub = {};
        mask_label = ind_l_tot;
        
        for k = 1:K-1,
            ind_temp = randperm(nb_labels - (k-1)*round(nb_labels/K),round(nb_labels/K));
            ind_l_sub{k} = mask_label(ind_temp);
            mask_label(ind_temp) = [];
        end
        ind_l_sub{K} = mask_label
        
        t = cputime;
        
        % color the selected snapshots as if they were unlabeled
        % dpERK
        snapshots_recolored1 = zeros(size(snapshots_orig3_1.images(:,:,3,:)));
        im_orig1 = double(snapshots_orig3_1.images(:,:,3,:));
        
        % ind
        snapshots_recolored2 = zeros(size(snapshots_orig3_1.images(:,:,2,:)));
        im_orig2 = double(snapshots_orig3_1.images(:,:,2,:));
        
        % rhomboid
        snapshots_recolored3 = zeros(size(snapshots_orig3_1.images(:,:,1,:)));
        im_orig3 = double(snapshots_orig3_1.images(:,:,1,:));
        
        % for each of the bins
        for k = 1:K
            indl_k = indl;
            
            indu_k = [indu,indl(ind_l_sub{k})];
            indl_k(ind_l_sub{k}) = [];
            
            l1 = length(indl_k);
            l2 = length(indu_k);
            
            % the semi-supervised learning problem is solved on each
            % combination of labeled and unlabeled data points
            [inv_u] = ssl_estimate( W, indl_k, indu_k);
            
            npixels = 512;
            
            e = cputime - t;
            fprintf([num2str(k),'/',num2str(K),' - Initialization - OK in ',num2str(e),' s \n']);
            t = cputime;
            
            for i = 1:npixels
                for j = 1:npixels
                    
                    Y1 = [double(squeeze(im_orig1(i,j,indl_k-min(indl_k)+1)))];
                    fu1 = inv_u*Y1;
                    
                    snapshots_recolored1(i,j,ind_l_sub{k}) = fu1(end-length(ind_l_sub{k})+1:end);
                    
                    Y2 = [double(squeeze(im_orig2(i,j,indl_k-min(indl_k)+1)))];
                    fu2 = inv_u*Y2;
                    
                    snapshots_recolored2(i,j,ind_l_sub{k}) = fu2(end-length(ind_l_sub{k})+1:end);
                    
                    Y3 = [double(squeeze(im_orig3(i,j,indl_k-min(indl_k)+1)))];
                    fu3 = inv_u*Y3;
                    
                    snapshots_recolored3(i,j,ind_l_sub{k}) = fu3(end-length(ind_l_sub{k})+1:end);
                end
            end
            
            e = cputime - t;
            fprintf(['Snapshots recolored in ',num2str(e),' s \n']);
            t = cputime;
            
        end
        
        % the recolored snapshots are compared to the original snapshots
        max_int = max(im_orig1(:));
        min_int = min(im_orig1(:));
        
        abs_error(u,r,1,3) = (1/length(im_orig1(:)))*sum(abs(snapshots_recolored1(:) - im_orig1(:)))/(max_int - min_int);
        
        max_int = max(im_orig2(:));
        min_int = min(im_orig2(:));
        
        abs_error(u,r,2,3) = (1/length(im_orig2(:)))*sum(abs(snapshots_recolored2(:) - im_orig2(:)))/(max_int - min_int);
        
        max_int = max(im_orig3(:));
        min_int = min(im_orig3(:));
        
        abs_error(u,r,3,3) = (1/length(im_orig3(:)))*sum(abs(snapshots_recolored3(:) - im_orig3(:)))/(max_int - min_int);
    end
end


%% K-fold cross validation and varying the number of unlabeled samples - data set 4
K = 2;

% varying the number of unlabeled samples
for u = 1:M
    fprintf([num2str(n_unlbds(u)),' \n']);
    
    % repetitions over the process of recoloring with K bins randomly drawn
    % and a fixed number of unlabeled samples
    for r = 1:nrep
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = find(all.movie_idx ~= 11 )';
        indl = find(all.movie_idx == 11 )';
        
        mask_unlbld = randperm(length(indu),nb_unlbld);
        indu_mask = sort(indu(mask_unlbld));
        mask_tot = [indu_mask,indl];
        nmov = length(indu_mask);
        
        W = AffinityFromDistance(V,10);
        
        % defining the K subsamples of the labeled data points
        nb_labels = length(snapshots_orig4.movie_idx);
        ind_l_tot = 1:nb_labels;
        
        ind_l_sub = {};
        mask_label = ind_l_tot;
        
        for k = 1:K-1
            ind_temp = randperm(nb_labels - (k-1)*round(nb_labels/K),round(nb_labels/K));
            ind_l_sub{k} = mask_label(ind_temp);
            mask_label(ind_temp) = [];
        end
        ind_l_sub{K} = mask_label;
        
        t = cputime;
        
        % color the selected snapshots as if they were unlabeled
        % Twist
        snapshots_recolored1 = zeros(size(snapshots_orig4_1.images(:,:,3,:)));
        im_orig1 = double(snapshots_orig4_1.images(:,:,3,:));
        
        % ind
        snapshots_recolored2 = zeros(size(snapshots_orig4_1.images(:,:,1,:)));
        im_orig2 = double(snapshots_orig4_1.images(:,:,1,:));
        
        % rhomboid
        snapshots_recolored3 = zeros(size(snapshots_orig4_1.images(:,:,2,:)));
        im_orig3 = double(snapshots_orig4_1.images(:,:,2,:));

        % for each of the bins
        for k = 1:K
            indl_k = indl;
            
            indu_k = [indu,indl(ind_l_sub{k})];
            indl_k(ind_l_sub{k}) = [];
            
            l1 = length(indl_k);
            l2 = length(indu_k);
            
            % the semi-supervised learning problem is solved on each
            % combination of labeled and unlabeled data points
            [inv_u] = ssl_estimate( W, indl_k, indu_k);
            
            npixels = 512;
            
            e = cputime - t;
            fprintf([num2str(k),'/',num2str(K),' - Initialization - OK in ',num2str(e),' s \n']);
            t = cputime;
            
            for i = 1:npixels
                for j = 1:npixels
                    
                    Y1 = [double(squeeze(im_orig1(i,j,indl_k-min(indl_k)+1)))];
                    fu1 = inv_u*Y1;
                    
                    snapshots_recolored1(i,j,ind_l_sub{k}) = fu1(end-length(ind_l_sub{k})+1:end);
                    
                    Y2 = [double(squeeze(im_orig2(i,j,indl_k-min(indl_k)+1)))];
                    fu2 = inv_u*Y2;
                    
                    snapshots_recolored2(i,j,ind_l_sub{k}) = fu2(end-length(ind_l_sub{k})+1:end);
                    
                    Y3 = [double(squeeze(im_orig3(i,j,indl_k-min(indl_k)+1)))];
                    fu3 = inv_u*Y3;
                    
                    snapshots_recolored3(i,j,ind_l_sub{k}) = fu3(end-length(ind_l_sub{k})+1:end);
                end
            end
            
            e = cputime - t;
            fprintf(['Snapshots recolored in ',num2str(e),' s \n']);
            t = cputime;
            
        end

        % the recolored snapshots are compared to the original snapshots
        max_int = max(im_orig1(:));
        min_int = min(im_orig1(:));
        
        abs_error(u,r,1,4) = (1/length(im_orig1(:)))*sum(abs(snapshots_recolored1(:) - im_orig1(:)))/(max_int - min_int);
        
        max_int = max(im_orig2(:));
        min_int = min(im_orig2(:));
        
        abs_error(u,r,2,4) = (1/length(im_orig2(:)))*sum(abs(snapshots_recolored2(:) - im_orig2(:)))/(max_int - min_int);
        
        max_int = max(im_orig3(:));
        min_int = min(im_orig3(:));
        
        abs_error(u,r,3,4) = (1/length(im_orig3(:)))*sum(abs(snapshots_recolored3(:) - im_orig3(:)))/(max_int - min_int);
    end
end

%%
d = date;
mkdir(d);
save([d,'/abs_error_all_datasets.mat'],'abs_error');


