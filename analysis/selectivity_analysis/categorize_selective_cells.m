function pools = categorize_selective_cells(selectivity_active, selectivity_passive, threshold, sig_cells)
    % Categorizes neurons into selectivity pools for both active and passive contexts
    %
    % Inputs:
    %   selectivity_active  - Vector of selectivity indices from active context
    %   selectivity_passive - Vector of selectivity indices from passive context
    %   threshold          - Numeric threshold for determining selectivity
    %                       (typically 0.1 or 0.2)
    %   sig_cells         - Indices of significantly modulated cells to analyze
    %
    % Output:
    %   pools - Structure with fields:
    %           .active.{left, right, nonsel}  - Cell indices for active context
    %           .passive.{left, right, nonsel} - Cell indices for passive context
    
    % Initialize output structure
    pools = struct();
    
    % Categorize cells in active context
    [pools.active.left, pools.active.right, pools.active.nonsel] = ...
        categorize_selectivity_context(selectivity_active, threshold, sig_cells);
    
    % Categorize cells in passive context
    [pools.passive.left, pools.passive.right, pools.passive.nonsel] = ...
        categorize_selectivity_context(selectivity_passive, threshold, sig_cells);
end