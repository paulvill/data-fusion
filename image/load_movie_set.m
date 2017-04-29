% LOAD_MOVIE_SET Loads a set of movies
%
% Usage
%    movies = load_movie_set(movie_opt);
%
% Input
%    movie_opt: An options structure specifying the location and format of the
%       movie files. The fields supported are:
%          - data_dir: The directory where the movies are found (required).
%          - image_name: The name of the image files, which is followed by a
%             number and extension to obtain the actual files (default
%             'image).
%          - image_ext: The extension of the image files (default 'tif').
%          - times_name: Name of the file containing times (default ...
%             'times.txt').
%          - angles_name: Name of the file containing angles (default 
%             'angles.txt').
%          - stack_name: The name of the stack (default '').
%          - nstack: The number of stacks (default 0).
%          - dim: The dimension of the images (default 2).
%       For more information on these options, see the documentation for
%       read_images.
%
% Output
%    movies: A structure containing the loaded movies, with fields:
%       - images: The actual images of the movies, arranged in a three-dimen-
%          sional array, where the third dimension is the image index.
%       - times: The time points of each image in the movies.
%       - angles: The angles through which the images have been rotated.
%       - movie_idx: The index of the movie each image belongs to.

function movies = load_movie_set(movie_opt)
    if ~isfield(movie_opt, 'data_dir')
        error('data_dir must be specified to load movie set');
    end

    movie_opt = fill_struct(movie_opt, ...
        'image_name', 'image', ...
        'image_ext', 'tif', ...
        'times_name', 'times.txt', ...
        'angles_name', 'angles.txt', ...
        'stack_name', '', ...
        'nstack', 0, ...
        'dim', 2, ...
        'npixels', [], ...
        'subset',[]);

    if numel(movie_opt.npixels) == 1
        movie_opt.npixels = movie_opt.npixels*ones(1, 2);
    end

    % Set up arrays to store images and associated data.
    images = {};
    times = {};
    angles = {};
    movie_idx = {};

    % Find movies (as subfolders or movie files).
    files = dir(movie_opt.data_dir);

    folders = files([files.isdir]==true);
    folders = {folders.name};

    folders = folders(~strncmp(folders, '.', 1));

    if numel(folders) == 0
        % No subdirectories are found so root must contain the images.
        folders = {''};
    end

    [images, times, angles, movie_idx, names] = ...
        load_image_movies(folders, movie_opt);
    
    % Reshape to get a flat structure for dataset.
    movies.images = cat(movie_opt.dim+1, images{:});
    movies.times = cat(1, times{:});
    movies.angles = cat(1, angles{:});
    movies.movie_idx = cat(1, movie_idx{:});
    movies.names = names(:);
end

function [images, times, angles, movie_idx, names] = load_image_movies( ...
    folders, movie_opt)
    movie_id = 1;
    for i = 1:numel(folders)
        image_path = fullfile(movie_opt.data_dir, folders{i});

        file_expr = [movie_opt.image_name '*.' movie_opt.image_ext];
        nimages = length(dir(fullfile(image_path, file_expr)));

        images{movie_id} = read_images(image_path, movie_opt.image_name, ...
            movie_opt.image_ext, movie_opt.stack_name, nimages, ...
            movie_opt.nstack, movie_opt.dim);

        if ~isempty(movie_opt.npixels)
            images{movie_id} = ...
                resize_images(images{movie_id}, movie_opt.npixels);
        end

        times_path = fullfile(image_path, movie_opt.times_name);
        if exist(times_path, 'file')
            times{movie_id} = dlmread(times_path);
        else
            times{movie_id} = zeros(nimages, 1);
        end

        angles_path = fullfile(image_path, movie_opt.angles_name);
        if exist(angles_path, 'file')
            angles{movie_id} = dlmread(angles_path);
        else
            angles{movie_id} = zeros(nimages, 1);
        end

        movie_idx{movie_id} = movie_id*ones(nimages, 1);

        names{movie_id} = folders{i};

        movie_id = movie_id+1;
    end
end
