function string = create_filename(first_string,type,N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,deltaN,estimated_noise_spectrum,NR_harmonics,estimate_transient)
    if isnan(deltaN)
        deltaN = "Inf";
    end
    string = first_string+"_"+type;
    values = [N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,NR_trials,RMS,flat_noise,deltaN,NR_harmonics];
    if strcmp(type,"WNLS")  
        values = [values,order,dof,estimate_transient,estimated_noise_spectrum];
    elseif strcmp(type,"TD")
        values = [values,l1];
    elseif strcmp(type,"FD")
        values = [values,l1,order,dof,estimate_transient];
    elseif strcmp(type,"FD_no_l1")
        values = [values,order,dof,estimate_transient];
    else
        error("Type not known")
    end
    for i = 1:numel(values)
        string = string + "_" + values(i);
    end
    string = string + ".mat";
end