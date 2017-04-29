% MKDIRP Recursive mkdir
%
% Usage
%    mkdirp(dir);
%
% Input
%    dir: A path to a directory to be created.
%
% Description
%    The function creates the directory and any parent directories that are
%    missing. Nothing is done for the directories already present. Mimics the
%    behavior of the -p flag to the GNU mkdir utility.

function mkdirp(dir)
    dir_parts = strsplit(dir, filesep);

    root = dir_parts{1};

    if ~exist(root, 'dir')
        mkdir(root);
    end

    if isempty(dir_parts{end})
        dir_parts = dir_parts(1:end-1);
    end

    for k = 2:numel(dir_parts)
        root = [root filesep dir_parts{k}];

        if ~exist(root, 'dir')
            mkdir(root);
        end
    end
end
