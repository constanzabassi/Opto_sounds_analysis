function combined_results = combine_dataset_results(all_results, total_cells_per_dataset)
% This function combines selectivity analysis results across multiple datasets and adjusts cell indices to be globally referenced across all datasets.
% Input Parameters:
% 
%   all_results: Cell array containing analysis results for each dataset
%   total_cells_per_dataset: Array containing the number of cells in each dataset
% 
% Processes two contexts: 'active' and 'passive'
% Handles three pool types: 'left', 'right', 'nonsel' (non-selective)
% Uses cumulative sum of cells to create global indices
%
% Output combined_results is resulting structure pooled across all datasets
% with fields 'active', 'passive', 'left','right','nonsel'

    % Combine results across datasets
    combined_results = struct();
    contexts = {'active', 'passive'};
    pool_types = {'left', 'right', 'nonsel'};
    
    % Calculate cumulative sum for global indexing
    cells_cumsum = [0; cumsum(total_cells_per_dataset)];
    
    for context_idx = 1:length(contexts)
        context = contexts{context_idx};
        for p_idx = 1:length(pool_types)
            pool_type = pool_types{p_idx};
            
            % Concatenate data across datasets
            left_mod = [];
            right_mod = [];
            max_mod = [];
            cell_indices = [];
            preferred_side = {};
            
            for dataset_idx = 1:length(all_results)
                dataset = all_results{dataset_idx};
                pool_data = dataset.(context).(pool_type);
                
                % Adjust cell indices to be global
                global_indices = pool_data.cell_indices + cells_cumsum(dataset_idx);
                
                left_mod = [left_mod; pool_data.left_mod(:)];
                right_mod = [right_mod; pool_data.right_mod(:)];
                max_mod = [max_mod; pool_data.max_mod(:)];
                cell_indices = [cell_indices; global_indices(:)];
                preferred_side = [preferred_side; pool_data.preferred_side(:)];
            end
            
            % Store combined results
            combined_results.(context).(pool_type) = struct(...
                'left_mod', left_mod, ...
                'right_mod', right_mod, ...
                'max_mod', max_mod, ...
                'preferred_side', {preferred_side}, ...
                'cell_indices', cell_indices);
            
            % Compute combined stats
            combined_results.(context).(pool_type).stats = compute_pool_stats(...
                combined_results.(context).(pool_type));
        end
    end
end