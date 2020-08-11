clear
close all
rng_index = 0;
rng(rng_index)

%% load configurations
study_case = "simple";
% study_case = "long_transient";
% study_case = "undermodeled";
for study_case = ["simple","long_transient","undermodeled"]
    clearvars -except rng_index study_case
    rng(rng_index)
    [keys,values] = load_config('base');
    for i = 1:numel(keys)
        evalc(keys{i} + "=[" + values{i}+"]");
    end

    [keys,values] = load_config(study_case);
    for i = 1:numel(keys)
        evalc(keys{i} + "=[" + values{i}+"]");
    end

    %% Simulate
    flat_noise = false;

    % System
    load(study_case + "/System_data")
    load(study_case + "/rho_actual")
    if ~exist("M","var")
        M = M_exact;
    end
    Fs = 1/Ts;
    qinv = tf([0,1],1,Ts,"variable","q^-1");

    % Frequency definition
    if even(N)
        k_max = N/2-1;
    else
        k_max = (N-1)/2;
    end

    % Signal
    Ntot = N*(P+P_SS+P_zero_input); %total number of samples
    t = (0:Ntot-1)*Ts;
    [r0,ExcitedHarm] = generate_MS(N,Fs,RMS,Fs/N,Fs/N*NR_harmonics,Fs/N);
    f_plot = ExcitedHarm*Fs/N;
    r0 = repmat(r0,P+P_SS,1);
    r0 = [zeros(N*P_zero_input,1);r0];
    d = randn(Ntot,NR_trials)*sigma_d;
    u0 = r0 + d;
    % output (noiseless)
    for i = NR_trials:-1:1
        y0(:,i) = lsim(G,u0(:,i),t);
    end

    % Noise filters
    if ~flat_noise
        load("Noise_filters")
    end

    % Noise
    ey_all = randn(Ntot,NR_trials);
    eu_all = randn(Ntot,NR_trials);
    ey_all = ey_all.*sigma_y;
    eu_all = eu_all.*sigma_u;

    for flat_noise = [true,false]
        if flat_noise
            vy_all = ey_all;
            vu_all = eu_all;
        else
            for i = NR_trials:-1:1
                vy_all(:,i) = lsim(Sy,ey_all(:,i),0:Ntot-1);
                vu_all(:,i) = lsim(Su,eu_all(:,i),0:Ntot-1);
            end
        end
        y_all = y0 + vy_all;
        u_all = u0 + vu_all;

        transient_analysis(y_all(:,1),P+P_SS+P_zero_input)

        % Transient removal
        r0 = r0((P_SS+P_zero_input)*N+1:end,:);
        u0 = u0((P_SS+P_zero_input)*N+1:end,:);
        y0 = y0((P_SS+P_zero_input)*N+1:end,:);
        vy_all = vy_all((P_SS+P_zero_input)*N+1:end,:);
        y_all = y_all((P_SS+P_zero_input)*N+1:end,:);
        vu_all = vu_all((P_SS+P_zero_input)*N+1:end,:);
        u_all = u_all((P_SS+P_zero_input)*N+1:end,:);
        t = t(1:end-(P_SS+P_zero_input)*N);
        Ntot = Ntot - (P_SS+P_zero_input)*N;

        %% Optimize
        % Controller
        for i = n_rho:-1:0
            beta_prime(i+1,1) = qinv^i; %better for some computations
            beta(i+1,1) = qinv^i / (1-qinv);
        end

        k_all = (0:N-1).';
        f_all = k_all*Fs/N;
        % frequence response for perfect filtering
        M_all = squeeze(freqresp(M,f_all,'Hz'));
        beta_all = squeeze(freqresp(beta,f_all,'Hz'));
        if n_rho > 0
            beta_all = beta_all.';
        end

        TD_method = TD_controller(P,M,F_all,beta,Phir_exact);
        FD_method = FD_controller(M_all,F_all,beta_all);
        FD_WNLS_method = FD_WNLS_controller(M_all,F_all,beta_all);

        for estimate_transient = [true,false]
            for type = ["TD", "FD", "WNLS", "FD_no_l1"]
                if strcmp(type,"WNLS") || strcmp(type,"FD_no_l1")
                    l1_array = 1;
                elseif strcmp(type,"FD")
                    l1_array = l1_array_FD;
                elseif strcmp(type,"TD")
                    l1_array = l1_array_TD;
                end
                for l1 = l1_array
                    string = create_filename(study_case,type,N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,deltaN,estimated_noise_spectrum,NR_harmonics,estimate_transient);
                    if isfile(study_case + "/data/" + string)
                        fprintf("Already exists\n")
                        continue
                    end
                    for trial_index = NR_trials+1-1:-1:1
                        y = y_all(:,trial_index);
                        u = u_all(:,trial_index);

                        % nonparametric estimate of the open loop plant
                        [G_est,varG] = nonparametric_estimate(r0,u,y,P,Ts,ExcitedHarm,order,dof,estimate_transient);

                        if strcmp(type,"WNLS")
                            % my method (WNLS)
                            [rho,~] = FD_WNLS_method.optimize(G_est,varG,ExcitedHarm,"LS",maxiter,precision);
                            rho_array_LS_init(:,trial_index) = rho;
                            [rho,~] = FD_WNLS_method.optimize(G_est,varG,ExcitedHarm,rho_actual,maxiter,precision);
                            rho_array_actual_init(:,trial_index) = rho;
                        elseif strcmp(type,"FD")
                            % my method (FD)
                            [rho,~] = FD_method.fast_optimize(G_est,ExcitedHarm,l1);
                            rho_array(:,trial_index) = rho;
                        elseif strcmp(type,"FD_no_l1")
                            [rho,~] = FD_method.fast_optimize_no_l1(G_est,ExcitedHarm);
                            rho_array(:,trial_index) = rho;
                        elseif strcmp(type,"TD")
                            % their method
                            [rho,~] = TD_method.fast_optimize(u,y,l1);
                            rho_array(:,trial_index) = rho;
                        else
                            error("Type doesn't exist")
                        end

                        fprintf("Done %d/%d\n",NR_trials-trial_index+1,NR_trials)
                    end

                    % Saving
                    if NR_trials >= 10
                        if strcmp(type,"WNLS")
                            save(study_case + "/data/" + string,...
                                'rho_array_LS_init','rho_array_actual_init','order','dof',...
                                'G','M','f_plot','N','P','P_SS','beta_prime','qinv','rng_index',...
                                'beta','ExcitedHarm')
                        elseif strcmp(type,"TD")
                            save(study_case + "/data/" + string,...
                                'rho_array','l1',...
                                'G','M','f_plot','N','P','P_SS','beta_prime','qinv','rng_index',...
                                'beta','ExcitedHarm')
                        elseif strcmp(type,"FD")
                            save(study_case + "/data/" + string,...
                                'rho_array','l1','order','dof',...
                                'G','M','f_plot','N','P','P_SS','beta_prime','qinv','rng_index',...
                                'beta','ExcitedHarm')
                        else
                            save(study_case + "/data/" + string,...
                                'rho_array','order','dof',...
                                'G','M','f_plot','N','P','P_SS','beta_prime','qinv','rng_index',...
                                'beta','ExcitedHarm')
                        end
                    end
                end
            end
        end
    end
end