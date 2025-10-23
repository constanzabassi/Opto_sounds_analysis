function T = struct2table_recursive(S, prefix, exclude_fields)
    if nargin < 2
        prefix = '';
    end
    if nargin < 3
        exclude_fields = {};
    end
    
    names = fieldnames(S);
    data = {};
    field_names = {};
    
    for i = 1:numel(names)
        field = names{i};
        fullfield = field;
        if ~isempty(prefix)
            fullfield = [prefix '.' field];
        end
        
        % Check if this field matches any exclusion
        is_excluded = any(cellfun(@(x) contains(fullfield, x), exclude_fields));
        if is_excluded
            continue; % skip this field
        end
        
        if isstruct(S.(field)) && ~contains(field,'p_va')
            % Recursively process sub-structure
            subT = struct2table_recursive(S.(field), fullfield, exclude_fields);
            field_names = [field_names; subT.FieldName];
            data = [data; subT.Value];
        elseif startsWith(field, 'p_val') 
            value = S.(field);
            % Handle different shapes of p_vals matrices
            [r, c] = size(value);

            if numel(value) == 1
                % Case ① single p-value
                field_names{end+1,1} = fullfield;
                data{end+1,1} = value;

            elseif r == 2
                % Case ② two rows: first = p-values, second = comparisons
                p_vals = value(1,:);
                comps  = value(2,:);
                for j = 1:numel(p_vals)
                    if j <= numel(comps)-1
                        name_suffix = sprintf('_%d_%d', comps{j});
                    else
                        name_suffix = sprintf('_%d', j);
                    end
                    field_names{end+1,1} = [fullfield name_suffix];
                    data{end+1,1} = p_vals{j};
                end

            elseif r == 1 && c > 2
                % Case ③ single row of multiple p-values
                for j = 1:c
                    field_names{end+1,1} = sprintf('%s_%d', fullfield, j);
                    data{end+1,1} = value{j};
                end

            elseif r > 2
                % Case ④ multi-row matrix with last row = comparison indices
                p_vals = value(1:end-1,:);
                comps  = value(end,:);
                num_sets = size(p_vals,1);
                for s = 1:num_sets
                    for j = 1:numel(comps)-1
                        name_suffix = sprintf('_%d_%d_%d', comps{j}, s);
                        field_names{end+1,1} = [fullfield name_suffix];
                        data{end+1,1} = p_vals{s,j};
                    end
                end

            else
                % Fallback case
                field_names{end+1,1} = fullfield;
                data{end+1,1} = value;
            end
        else
            % Store field name and value
            field_names{end+1,1} = fullfield;
            data{end+1,1} = S.(field);
        end
    end
    
    T = table(field_names, data, 'VariableNames', {'FieldName', 'Value'});
end
