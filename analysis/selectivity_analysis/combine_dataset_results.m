function combined_results = combine_dataset_results(all_results, total_cells_per_dataset, all_celltypes)
    % Initialize result structures
    combined_results = struct();
    contexts = {'active', 'passive','both'};
    pool_types = {'left', 'right', 'nonsel'};
    
    % Calculate cumulative sum for global indexing
    cells_cumsum = [0; cumsum(total_cells_per_dataset)];
    
    % Process data for each context and pool type
    [combined_results, all_pool_indices] = process_pool_data(all_results, cells_cumsum, contexts, pool_types);
    
    % Add cell type information if provided
    if nargin > 2 && ~isempty(all_celltypes)
        combined_results.celltypes = categorize_celltypes_by_pool(all_pool_indices, all_celltypes, cells_cumsum);
    end
end

function [combined_results, all_pool_indices] = process_pool_data(all_results, cells_cumsum, contexts, pool_types)
    % Process and combine pool data across datasets
    %
    % Inputs:
    %   all_results    - Cell array of results from each dataset
    %   cells_cumsum   - Cumulative sum of cells for global indexing
    %   contexts       - Cell array of context names
    %   pool_types     - Cell array of pool type names
    %
    % Outputs:
    %   combined_results  - Structure with combined data across datasets
    %   all_pool_indices - Structure containing global indices for each pool
    
    % Initialize structures
    combined_results = struct();
    all_pool_indices = struct();
    
    % Process each context and pool type
    for ctx_idx = 1:length(contexts)
        context = contexts{ctx_idx};
        for p_idx = 1:length(pool_types)
            pool_type = pool_types{p_idx};
            
            % Initialize arrays for current pool
            pool_data = struct(...
                'left_mod', [], ...
                'right_mod', [], ...
                'max_mod', [], ...
                'preferred_side', {{}}, ...
                'cell_indices', []);
            indices = [];
            
            % Concatenate data across datasets
            for d_idx = 1:length(all_results)
                dataset = all_results{d_idx};
                curr_pool = dataset.(context).(pool_type);
                
                % Adjust indices for global reference
                global_indices = curr_pool.cell_indices + cells_cumsum(d_idx);
                
                % Accumulate data
                pool_data = accumulate_pool_data(pool_data, curr_pool, global_indices);
                indices = [indices; global_indices(:)];
            end
            
            % Store combined results
            combined_results.(context).(pool_type) = create_pool_struct(pool_data);
            all_pool_indices.(context).(pool_type) = indices;
        end
    end
end

function pool_data = accumulate_pool_data(pool_data, curr_pool, global_indices)
    pool_data.left_mod = [pool_data.left_mod; curr_pool.left_mod(:)];
    pool_data.right_mod = [pool_data.right_mod; curr_pool.right_mod(:)];
    pool_data.max_mod = [pool_data.max_mod; curr_pool.max_mod(:)];
    pool_data.preferred_side = [pool_data.preferred_side; curr_pool.preferred_side(:)];
    pool_data.cell_indices = [pool_data.cell_indices; global_indices(:)];  % Added this line
end

function pool_struct = create_pool_struct(pool_data)
    pool_struct = struct(...
        'left_mod', pool_data.left_mod, ...
        'right_mod', pool_data.right_mod, ...
        'max_mod', pool_data.max_mod, ...
        'preferred_side', {pool_data.preferred_side}, ...
        'cell_indices', pool_data.cell_indices);
    
    % Add statistics
    pool_struct.stats = compute_pool_stats(pool_struct);
end

function celltype_pools = categorize_celltypes_by_pool(all_pool_indices, all_celltypes, cells_cumsum)
    % Initialize cell type tracking structure
    celltype_pools = struct();
    celltypes = {'pyr', 'som', 'pv'};
    contexts = {'active', 'passive','both'};
    pool_types = {'left', 'right', 'nonsel'};
    
    for ctx_idx = 1:length(contexts)
        context = contexts{ctx_idx};
        for pool_idx = 1:length(pool_types)
            pool = pool_types{pool_idx};
            
            % Get global indices for current pool
            pool_indices = all_pool_indices.(context).(pool);
            
            % Categorize by cell type
            for type_idx = 1:length(celltypes)
                ctype = celltypes{type_idx};
                celltype_pools.(ctype).(context).(pool) = ...
                    find_celltype_in_pool(pool_indices, all_celltypes, ctype, cells_cumsum);
            end
        end
    end
end

function type_indices = find_celltype_in_pool(pool_indices, all_celltypes, celltype, cells_cumsum)
    type_indices = [];
    field_name = [celltype '_cells'];
    
    for dataset_idx = 1:length(all_celltypes)
        % Get cells of this type in current dataset
        dataset_cells = all_celltypes{1,dataset_idx}.(field_name);
        
        % Adjust indices for global reference
        global_cells = dataset_cells + cells_cumsum(dataset_idx);
        
        % Find intersection with pool indices
        type_indices = [type_indices; intersect(pool_indices, global_cells)];
    end
end