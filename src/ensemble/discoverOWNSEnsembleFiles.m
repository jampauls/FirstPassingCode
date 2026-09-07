function fileList = discoverOWNSEnsembleFiles(directoryPath)
%DISCOVEROWNSENSEMBLEFILES List all MAT-files in a directory.
%
%   fileList = discoverOWNSEnsembleFiles(directoryPath)
%
% All MAT-files in directoryPath are assumed to be OWNS solution files for
% one (omega,beta) mode each, regardless of filename. Files are returned
% in sorted order for reproducibility.
%
% MATLAB version: R2020b

if nargin < 1 || isempty(directoryPath)
    error('discoverOWNSEnsembleFiles:MissingDirectory', ...
        'directoryPath is required.');
end

if ~exist(directoryPath, 'dir')
    error('discoverOWNSEnsembleFiles:DirectoryNotFound', ...
        'Directory not found: %s', directoryPath);
end

listing = dir(fullfile(directoryPath, '*.mat'));
listing = listing(~[listing.isdir]);

[~, order] = sort({listing.name});
listing = listing(order);

if isempty(listing)
    error('discoverOWNSEnsembleFiles:NoFilesFound', ...
        'No MAT-files were found in: %s', directoryPath);
end

fileList = cell(numel(listing), 1);

for k = 1:numel(listing)
    fileList{k} = fullfile(listing(k).folder, listing(k).name);
end

end
