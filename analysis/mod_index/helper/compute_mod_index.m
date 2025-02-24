function mod_index = compute_mod_index(avg1, avg2, mod_type)
%depending on mode computes modulatio index
    switch mod_type
        case 'ctrl'
            mod_index = compute_mod_index_ctrl(avg1, avg2);
        case 'influence'
            mod_index = compute_mod_index_influence(avg1, avg2);
        case 'prepost'
            mod_index = compute_mod_index_prepost(avg1, avg2); %post,pre
        case 'prepost_sound'
            mod_index = compute_mod_index_prepost(avg1, avg2); %post,pre 
        otherwise
            error('Invalid mod_type. Choose from ''ctrl'', ''influence'', or ''prepost''');
    end
end