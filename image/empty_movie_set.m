% EMPTY_MOVIE_SET Construct an empty movie set
%
% Usage
%    movie_set = empty_movie_set(npixels, channels);
%
% Input
%    npixels: The width (and height) of the images in the (square) movies
%       (default 100).
%    channels: The number of channels in the movies (default 1).
%
% Output
%    movie_set: An empty set of movies.

function movie_set = empty_movie_set(npixels, channels)
    if nargin < 1 || isempty(npixels)
        npixels = 100;
    end

    if nargin < 2 || isempty(channels)
        channels = 1;
    end

    movie_set = struct();

    if channels == 1
        movie_set.images = zeros([npixels*ones(1, 2) 0]);
    else
        movie_set.images = zeros([npixels*ones(1, 2) channels 0]);
    end

    movie_set.times = zeros(0, 1);
    movie_set.angles = zeros(0, 1);
    movie_set.movie_idx = zeros(0, 1);
    movie_set.names = {};
end
