% CAT_MOVIE_SETS Concatenate movie sets
%
% Usage
%    movie_set = cat_movie_sets(movie_set1, movie_set2);
%
% Input
%    movie_set1, movie_set2: Two movie sets to be concatenated.
%
% Output
%    movie_set: The concatenation of the two movie sets.

function movie_set = cat_movie_sets(movie_set1, movie_set2)
    ndims1 = ndims(movie_set1.images);
    ndims2 = ndims(movie_set2.images);

    if ndims1 ~= ndims2 ...
        || (ndims1 == ndims2 && ndims1 == 4 ...
            && size(movie_set1.images, 3) ~= size(movie_set2.images, 3))
        error('number of channels must be the same');
    end

    if ndims1 == 3
        movie_set1.images = permute(movie_set1.images, [1 2 4 3]);
        movie_set2.images = permute(movie_set2.images, [1 2 4 3]);
    end

    movie_set = movie_set1;
    movie_set.images = cat(4, movie_set.images, movie_set2.images);
    movie_set.times = cat(1, movie_set.times, movie_set2.times);
    movie_set.angles = cat(1, movie_set.angles, movie_set2.angles);
    movie_set.movie_idx = cat(1, movie_set.movie_idx, ...
        max([0; movie_set.movie_idx])+movie_set2.movie_idx);
    movie_set.names = cat(1, movie_set.names, movie_set2.names);

    if ndims1 == 3
        movie_set.images = permute(movie_set.images, [1 2 4 3]);
    end
end

