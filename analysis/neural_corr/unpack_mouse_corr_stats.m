function unpacked_index_matrix = unpack_mouse_corr_stats(mouse_corr_stats, mode)
    % Extract velocity types
    field_types = fieldnames(mouse_corr_stats);
        
    %Get total contexts in the structure
    total_contexts = size(mouse_corr_stats.(field_types{1}),2);

    corr_data_all = [];
    unpacked_index_matrix =[];
    
    % Loop through velocity types
    for field_idx = 1:length(field_types)
        field_type = field_types{field_idx};
        
        % Loop through contexts
        for context_idx = 1:total_contexts
            
            
                for direction = 1:2 %left and right
                    data = []; % Collect correlation data across mice
                    % Loop through mice
                    for mouse_idx = 1:size(mouse_corr_stats.(field_type), 1)
                        data = [data; mouse_corr_stats.(field_type){mouse_idx, context_idx, direction}(:)]; % Use first or second index based on requirement
                    end
    
                %Save across mice
                corr_data_all(field_idx,context_idx,:,direction) = data;


                end
                
                % Get max or mean across left and right
                if strcmp(mode,'max')
                    unpacked_index_matrix(field_idx,:,context_idx) = max(abs(corr_data_all(field_idx,context_idx,:,1)),abs(corr_data_all(field_idx,context_idx,:,2)));
                elseif strcmp(mode,'mean')
                    unpacked_index_matrix(field_idx,:,context_idx) = mean([squeeze([corr_data_all(field_idx,context_idx,:,1)]),squeeze([corr_data_all(field_idx,context_idx,:,2)])],2);

                end

                unpacked_index_matrix = squeeze(unpacked_index_matrix); %if only one field it should squeeze it out

        end
    end
