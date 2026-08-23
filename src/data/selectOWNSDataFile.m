function dataFile = selectOWNSDataFile(ownsCfg)
%SELECTOWNSDATAFILE Select an OWNS solution file using UIGETFILE.
%
%   dataFile = selectOWNSDataFile(ownsCfg)
%
% Supported configuration fields:
%
%   ownsCfg.userSelectDirectory
%       Initial directory displayed by uigetfile.
%
%   ownsCfg.userSelectTitle
%       File-selection dialog title.
%
%   ownsCfg.userSelectDefaultFile
%       Optional default filename.
%
%   ownsCfg.userSelectVerbose
%       Print the selected full path.
%
% Output:
%   dataFile
%       Full path to the selected MAT-file.
%
% An error is raised if the user cancels the dialog.
%
% MATLAB version: R2020b

if nargin < 1 || isempty(ownsCfg)
    ownsCfg = struct();
end

initialDirectory = getOption( ...
    ownsCfg, ...
    'userSelectDirectory', ...
    pwd);

dialogTitle = getOption( ...
    ownsCfg, ...
    'userSelectTitle', ...
    'Select an OWNS solution MAT-file');

defaultFile = getOption( ...
    ownsCfg, ...
    'userSelectDefaultFile', ...
    '');

verbose = getOption( ...
    ownsCfg, ...
    'userSelectVerbose', ...
    true);

% Convert string values to character arrays for MATLAB R2020b
% compatibility.
if isstring(initialDirectory)
    initialDirectory = char(initialDirectory);
end

if isstring(dialogTitle)
    dialogTitle = char(dialogTitle);
end

if isstring(defaultFile)
    defaultFile = char(defaultFile);
end

if isempty(initialDirectory)
    initialDirectory = pwd;
end

% Resolve a relative directory against the current working directory.
if ~isAbsolutePath(initialDirectory)
    initialDirectory = fullfile(pwd, initialDirectory);
end

if ~exist(initialDirectory, 'dir')
    warning('selectOWNSDataFile:InitialDirectoryNotFound', ...
        ['The configured OWNS selection directory does not exist:\n', ...
         '  %s\n', ...
         'The current working directory will be used instead.'], ...
        initialDirectory);

    initialDirectory = pwd;
end

% The third argument to uigetfile can be a full default path. If no
% default filename was supplied, use "*.mat" so the dialog opens in the
% requested directory and initially displays MAT-files.
if isempty(defaultFile)
    defaultPath = fullfile(initialDirectory, '*.mat');
else
    defaultPath = fullfile(initialDirectory, defaultFile);
end

filterSpecification = { ...
    '*.mat', 'MAT-files (*.mat)'; ...
    '*.*',   'All files (*.*)'};

[selectedFile, selectedPath] = uigetfile( ...
    filterSpecification, ...
    dialogTitle, ...
    defaultPath, ...
    'MultiSelect', 'off');

% uigetfile returns numeric zero when the user cancels.
if isequal(selectedFile, 0) || isequal(selectedPath, 0)
    error('selectOWNSDataFile:SelectionCancelled', ...
        'OWNS file selection was cancelled by the user.');
end

dataFile = fullfile(selectedPath, selectedFile);

if ~isfile(dataFile)
    error('selectOWNSDataFile:SelectedFileNotFound', ...
        'The selected file does not exist: %s', dataFile);
end

if verbose
    fprintf('Selected OWNS solution file:\n');
    fprintf('  %s\n', dataFile);
end

end

% =========================================================================
function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end

% =========================================================================
function result = isAbsolutePath(pathName)
%ISABSOLUTEPATH Test whether a path is absolute on Windows or UNIX.

if isempty(pathName)
    result = false;
    return;
end

if ispc
    % Match paths such as:
    %
    %   C:\data
    %   C:/data
    %   \\server\share
    result = ...
        (~isempty(regexp(pathName, '^[A-Za-z]:[\\/]', 'once'))) ...
        || startsWith(pathName, '\\');
else
    result = startsWith(pathName, '/');
end

end
