function [x, B, H, eThresh, meta] = loadProblemData(dataFile, cfg)
%LOADPROBLEMDATA Load precomputed first-transition input data.

if ~isfile(dataFile)
    error('loadProblemData:FileNotFound', ...
        'Data file not found: %s', dataFile);
end

data = load(dataFile);

required = {'x', 'B', 'H', 'eThresh'};

for j = 1:numel(required)
    if ~isfield(data, required{j})
        error('loadProblemData:MissingVariable', ...
            'The data file does not contain "%s".', required{j});
    end
end

x = data.x;
B = data.B;
H = data.H;
eThresh = data.eThresh;

if isfield(data, 'meta')
    meta = data.meta;
else
    meta = struct();
end

meta.sourceFile = dataFile;
meta.cfgAtLoad = cfg;

end