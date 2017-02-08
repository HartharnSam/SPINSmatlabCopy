function fields = find_fields(ind)
% FIND_FIELDS     Find the SPINS fields for a given index (ind)
%                 in the working directory
%
%  Usage:
%    fields = find_fields()
%    fields = find_fields(5)
%
%  Inputs:
%    - one optional argument. A integer of the output index to use
%
%  Outputs:
%    'fields' - cell array of all fields for given output index
%
%  David Deepwell, 2017

% use the index to find the fields
    if exist('ind','var')
        files = dir(['*.',num2str(ind)]);
    else
        files = dir('*.0');
    end
    nfiles = length(files);
    table = struct2dataset(files);

    % find all field names
    fields = cell(1,nfiles);
    for ii = 1:nfiles
        filename = table{ii,1};
        [field, ~] = strtok(filename, '.');
        fields{ii} = field;
    end
end
