clear all
% close all

%% read in data

% red channel contains signal of interest
channel = 1;

% subsample images 
npixels = 100;

% timestep for movies
dt = 0.5;

movie_path = 'data/movies';

% live imaging movies
movies = {'bomyi_emb01_gast01.avi';...
    'bomyi_emb02_gast02.avi';...
    '14_0623_emb01_hisRFP_gastrulation.avi';...
    '14_0623_emb02_hisRFP_gastrulation.avi';...
    '14_0624_emb01_hisRFP_gastrulation.avi';...
    '14_0624_emb02_hisRFP_gastrulation.avi'; ...
    '0709_emb02_cell_gast.avi'};

movies = cellfun(@(f)(fullfile(movie_path, f)), movies, ...
	'UniformOutput', false);

nmovies = length(movies);

% rotation angle for each movie so that the dorsal side is at the bottom
theta = [-95;
    -55;
    -80;
    -80;
    -95;
    -120;
    -95];


% start and end frames for each movie so that they are (approximately)
% temporally aligned
movie_start = [20 25 15 10 15 5 110];
  movie_end = [55 60 48 42 42 40 140];

% function to make subplots
make_subplot = @(i) subplot(2, 4, i);

% store images and times (each movie is stored in one cell of the cell
% array)
images = cell(nmovies, 1);
time = cell(nmovies, 1);
nimages = zeros(nmovies, 1);

% number of PCA modes to use
nmodes = 5;

%% read in images
figure;
for i=1:nmovies
    
    % read in movie
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
    
    % extract selected frames
    images_tmp = images_tmp(:, :, :, movie_start(i):movie_end(i));
    
    % calculate developmental times
    time_tmp = time_tmp(movie_start(i):movie_end(i)) - time_tmp(movie_start(i));
    time{i} = time_tmp * dt;
    
    images_tmp = images_tmp(:, :, channel, :);
    images_tmp = squeeze(images_tmp);
    
    nimages(i) = length(time{i});
    
    % normalize images from movie (normalize signal intensity and crop)
    for j=1:nimages(i)
        images_tmp(:, :, j) = image_normalize(images_tmp(:, :, j), ...
            theta(i), true);
    end
    
    images{i} = images_tmp;
    
    % plot representative image 
    make_subplot(i);
    imagesc(images{i}(:,:, round(nimages(i)/2)));
    colormap('gray');
    set(gca,'XTick',[],'YTick',[]);
    axis square
end

%% unwrap imaging data into a 2D matrix for PCA

PCA_data = cell(nmovies, 1);
for i=1:nmovies
    PCA_data{i} = create_PCA_data(images{i});
end

%% calculate shift times to better temporally align movis

shift_times = estimate_shift_times(PCA_data, time, nmodes);

time_adjust = time;
for i=1:nmovies
    time_adjust{i} = time_adjust{i} - shift_times(i);
end

%% fit linear model (time as a function of PCA coefficients)

train_data = vertcat(PCA_data{:});
train_times = vertcat(time_adjust{:});

model = train_PCA_reg_model(train_data, train_times, nmodes);

% save data
save('PCA_train_data.mat', '-struct', 'model');

%% plot fit vs. estimated times

for i=1:nmovies
    [images_tmp, time_tmp] = read_video(movies{i}, npixels);
    images_tmp = images_tmp(:, :, :, movie_start(i):movie_end(i));
    
    time_est = zeros(size(time_adjust{i}));
    time_est_min = zeros(size(time_adjust{i}));
    time_est_max = zeros(size(time_adjust{i}));
    
    for j=1:nimages(i)
        images_j = imrotate(images_tmp(:,:,channel,j), theta(i), 'crop');
        [time_est(j),time_est_min(j),time_est_max(j)] = ...
            predict_time_image(model, images_j);
    end
    
    figure;
    subplot(1,3,1)
    plot(time_adjust{i}, time_est,'.')
    hold on
    plot([0 15], [0 15],'-r')

    subplot(1,3,2)
    plot(abs(time_adjust{i} - time_est),'.')

    subplot(1,3,3)
     plot(time_est_max - time_est_min)

end

