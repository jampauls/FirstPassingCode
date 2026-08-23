function Bsyn = convertFactorOrientation(Bin, factorOrientation)
%CONVERTFACTORORIENTATION Convert factor data to synthesis orientation.
%
% Internal convention:
%
%   q = Bsyn*w
%   C = Bsyn*Bsyn'
%
% Supported input conventions:
%
%   'synthesis':
%       C = B*B'
%
%   'row':
%       C = B_row'*B_row
%       Bsyn = B_row'
%
% Bin may be a cell array or a single numeric matrix.

if nargin < 2 || isempty(factorOrientation)
    factorOrientation = 'synthesis';
end

if isstring(factorOrientation)
    factorOrientation = char(factorOrientation);
end

switch lower(factorOrientation)

    case 'synthesis'
        Bsyn = Bin;

    case 'row'
        if iscell(Bin)
            Bsyn = cell(size(Bin));

            for n = 1:numel(Bin)
                Bsyn{n} = Bin{n}';
            end
        elseif isnumeric(Bin)
            Bsyn = Bin';
        else
            error('convertFactorOrientation:InvalidInput', ...
                'Bin must be a cell array or numeric matrix.');
        end

    otherwise
        error('convertFactorOrientation:UnknownOrientation', ...
            'Unknown factor orientation "%s".', factorOrientation);
end

end
