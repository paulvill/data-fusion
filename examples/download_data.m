% DOWNLOAD_DATA Download images for example scripts
%
% Usage
%     download_data()
%
% Description
%     Downloads the images used for the example scripts experimental_dataset.m
%     and experimental_dataset_cross_validation.m. The images are saved in the
%     data subdirectory. The amount of space needed is 643 MB.

function download_data()
    url = 'https://github.com/paulvill/data-fusion-images/archive/master.zip';

    location = fullfile('data', 'data-fusion-images.zip');

    mkdirp('data');

    fprintf('Downloading data into %s...', location);
    urlwrite(url, location);
    fprintf('OK\n');

    fprintf('Extracting archive...');
    filenames = unzip(location, 'data');
    delete(location);

    for k = 1:numel(filenames)
        source = filenames{k};

        [source_path, source_file, source_ext] = fileparts(source);
        [~, source_path] = fileparts(source_path);

        dest = fullfile('data', source_path, [source_file source_ext]);

        mkdirp(fileparts(dest));

        if ~strcmp(source, dest)
            movefile(source, dest);
        end
    end

    if isoctave()
        confirm_recursive_rmdir(false, 'local');
    end
    rmdir(fullfile('data', 'data-fusion-images-master'), 's');
    fprintf('OK\n');
end
