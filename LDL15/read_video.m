function [images, times] = read_video(avi_file, npixels)

if nargin < 2
    npixels = [];
end

if numel(npixels) == 1
    npixels = npixels*ones(1, 2);
end

vidobj = VideoReader(avi_file);

num_images = vidobj.NumberOfFrames;

if isempty(npixels)
    npixels = [vidobj.Height vidobj.Width];
end

nchannels = 3;
images = zeros(npixels(1), npixels(2), nchannels, num_images, 'uint8');

for i = 1:num_images
    A = read(vidobj, i);
    A = imresize(A, npixels);
    images(:, :, :, i) = A;
end

times = (1:num_images)';

