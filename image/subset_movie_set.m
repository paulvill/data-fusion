% SUBSET_MOVIE_SET Extract a subset of the movie images
%
% Usage
%    movies = subset_movie_set(movies, mask);
%
% Input
%    movies: The movies dataset from which the subset is to be extracted.
%    mask: The mask specifiying which images to extract.
%
% Output
%    movies: The set of movies corresponding to the subset of images specified
%       in mask.

function movies = subset_movie_set(movies, mask)
    dim = ndims(movies.images)-1;

    % Depending on dimension, subset last index.
    if dim == 2
        movies.images = movies.images(:,:,mask);
    elseif dim == 3
        movies.images = movies.images(:,:,:,mask);
    else
        error('dimension must be 2 or 3');
    end

    movies.times = movies.times(mask);
    movies.angles = movies.angles(mask);
    movies.movie_idx = movies.movie_idx(mask);

    % Reindex the movie_idx field to be contiguous integers.
    [~,~,movies.movie_idx] = unique(movies.movie_idx);
end

