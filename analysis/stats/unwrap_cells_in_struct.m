function S = unwrap_cells_in_struct(S)
% UNWRAP_CELLS_IN_STRUCT Recursively unwraps single-cell contents of structs.
% 
% Example:
%   S.pos_percents{1,1} = struct('mean', 0.5, 'sd', 0.1);
%   S = unwrap_cells_in_struct(S);
%   % Now:
%   % S.pos_percents.mean = 0.5
%   % S.pos_percents.sd   = 0.1

    fields = fieldnames(S);
    for i = 1:numel(fields)
        f = fields{i};
        val = S.(f);

        % 🔹 If value is a cell with a single struct inside, unwrap it
        if iscell(val) && numel(val) == 1 && isstruct(val{1})
            % Replace the cell with the struct
            S.(f) = val{1};
            % Recurse into that struct in case it has cells inside too
            S.(f) = unwrap_cells_in_struct(S.(f));

        % 🔹 If value is a cell array of structs (e.g., {1xN})
        elseif iscell(val) && all(cellfun(@(x) isstruct(x), val))
            % Merge cell structs into indexed subfields
            tmp = struct();
            for j = 1:numel(val)
                subname = sprintf('entry%d', j);
                tmp.(subname) = unwrap_cells_in_struct(val{j});
            end
            S.(f) = tmp;

        % 🔹 If it’s a struct itself, go deeper recursively
        elseif isstruct(val)
            S.(f) = unwrap_cells_in_struct(val);
        end
    end
end
