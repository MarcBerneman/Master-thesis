clear
close all
save_figs = false;

study_case = "undermodeled";

%% load configurations
[keys,values] = load_config('base');
for i = 1:numel(keys)
    evalc(keys{i} + "=[" + values{i}+"]");
end

[keys,values] = load_config(study_case);
for i = 1:numel(keys)
    evalc(keys{i} + "=[" + values{i}+"]");
end
clearvars keys values

flat_noise = false;

%% create maps
titles_map = containers.Map;
array_map = containers.Map;
J_map = containers.Map;
Jmr_map = containers.Map;
mean_J_map = containers.Map;
mean_Jmr_map = containers.Map;
RMS_map = containers.Map;
CL_map = containers.Map;
l1_map = containers.Map;

%% load data
for i = numel(l1_array_TD):-1:1
    l1 = l1_array_TD(i);
    string = create_filename(study_case,"TD",N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,NaN,estimated_noise_spectrum,NR_harmonics,estimate_transient);
    load(study_case + "/data/" + string,"rho_array")
    array_map("TD_" + l1) = rho_array;
    titles_map("TD_" + l1) = char("TD (l_1 = " + l1 + ")");
    l1_map("TD_" + l1) = l1;
end

for estimate_transient = [true,false]
    for i = numel(l1_array_FD):-1:1
        l1 = l1_array_FD(i);
        string = create_filename(study_case,"FD",N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,NaN,estimated_noise_spectrum,NR_harmonics,estimate_transient);
        load(study_case + "/data/" + string,"rho_array")
        array_map("FD_" + estimate_transient + "_" + l1) = rho_array;
        if estimate_transient
            titles_map("FD_" + estimate_transient + "_" + l1) = char("FD (l_1 = " + l1 + ", transient)");
        else
            titles_map("FD_" + estimate_transient + "_" + l1) = char("FD (l_1 = " + l1 + ", no transient)");
        end
        l1_map("FD_" + estimate_transient + "_" + l1) = l1;
    end

    string = create_filename(study_case,"FD_no_l1",N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,NaN,estimated_noise_spectrum,NR_harmonics,estimate_transient);
    load(study_case + "/data/" + string,"rho_array")
    array_map("FD_" + estimate_transient + "_no_l1") = rho_array;
    if estimate_transient
        titles_map("FD_" + estimate_transient + "_no_l1") = char("FD (transient)");
    else
        titles_map("FD_" + estimate_transient + "_no_l1") = char("FD (no transient)");
    end
    l1_map("FD_" + estimate_transient + "_no_l1") = N-1;

    string = create_filename(study_case,"WNLS",N,P,P_SS,P_zero_input,sigma_u,sigma_y,sigma_d,l1,order,dof,NR_trials,RMS,flat_noise,NaN,estimated_noise_spectrum,NR_harmonics,estimate_transient);
    load(study_case + "/data/" + string,"rho_array_LS_init")
    load(study_case + "/data/" + string,"rho_array_actual_init")
    array_map("WNLS_" + estimate_transient + "_actual") = rho_array_actual_init;
    array_map("WNLS_" + estimate_transient + "_LS") = rho_array_LS_init;
    if estimate_transient
        titles_map("WNLS_" + estimate_transient + "_actual") = char("WNLS (optimal init, transient)");
        titles_map("WNLS_" + estimate_transient + "_LS") = char("WNLS (LS init, transient)");
    else
        titles_map("WNLS_" + estimate_transient + "_actual") = char("WNLS (optimal init, no transient)");
        titles_map("WNLS_" + estimate_transient + "_LS") = char("WNLS (LS init, no transient)");
    end
end
load(study_case + "/rho_actual","rho_actual")
array_map("actual") = rho_actual;
titles_map("actual") = char("Optimal");

load(study_case + "/data/" + string,"M","G","beta","ExcitedHarm")

if even(N)
    k_max = N/2-1;
else
    k_max = (N-1)/2;
end
Fs = 1/M.Ts;
f = (0:k_max)*Fs/N;
f_exc = f(ExcitedHarm+1);
Gk = squeeze(freqresp(G,f_exc,'Hz'));
betak = squeeze(freqresp(beta,f_exc,'Hz'));
Mk = squeeze(freqresp(M,f_exc,'Hz'));

keystrings = keys(array_map);
for i = 1:numel(keystrings)
    key = keystrings{i};
    [Jmr,J] = fast_calc_cost(Mk,Gk,betak,array_map(key));
    Jmr_map(key) = Jmr;
    J_map(key) = J;
    mean_Jmr_map(key) = mean(Jmr);
    mean_J_map(key) = mean(J);
    RMS_map(key) = fast_calc_RMS_diff(Mk,Gk,betak,array_map(key));
    CL_map(key) = fast_calc_CL(Gk,betak,array_map(key));
end

%% figures
tit1 = "Mean cost function at excited frequencies (dB) | " + NR_trials + " trials | TD";
if flat_noise
    tit2 = "White noise";
else
    tit2 = "Coloured noise";
end

figure('Units','Normalized','Position',[0.1728    0.4688    0.4876    0.3919])
hold on
plot(l1_array_TD,0.5*db(getValues(mean_Jmr_map,"TD_" + l1_array_TD)),'b-','Linewidth',4,'Markersize',15,'Displayname','J_{mr}')
plot(l1_array_TD,0.5*db(getValues(mean_J_map,"TD_" + l1_array_TD)),'r:','Linewidth',4,'Markersize',15,'Displayname','J')
xlabel('l_1')
legend("Location","Best")
title([tit1,tit2])
set(gca,'xscale','log')
set(gca,'YMinorGrid','off')
set(gca,'XMinorGrid','off')
xlim([1,max(l1_array_TD)])
xticklabels([1,10,100])
plot_options(gca)
if strcmp(study_case,"simple") && flat_noise && save_figs
    print(gcf,"report_figures/figures/mean_cost_function_TD_simple_flat","-depsc")
end

estimate_transient = true;

tit1 = "Mean cost function at excited frequencies (dB) | " + NR_trials + " trials | FD";
tit2 = "Transient";
if ~estimate_transient
    tit2 = tit2 + " not";
end
tit2 = tit2 + " suppressed | ";
if flat_noise
    tit2 = tit2 + "white noise";
else
    tit2 = tit2 + "coloured noise";
end

figure('Units','Normalized','Position',[0.1728    0.4688    0.4876    0.3919])
hold on
plot(l1_array_FD,0.5*db(getValues(mean_Jmr_map,"FD_" + estimate_transient + "_" + l1_array_FD)),'b-','Linewidth',4,'Markersize',15,'Displayname','J_{mr}')
plot(l1_array_FD,0.5*db(getValues(mean_J_map,"FD_" + estimate_transient + "_" + + l1_array_FD)),'r:','Linewidth',4,'Markersize',15,'Displayname','J')
xlabel('l_1')
legend("Location","Best")
title({tit1,tit2})
set(gca,'xscale','log')
set(gca,'YMinorGrid','off')
set(gca,'XMinorGrid','off')
xlim([1,max(l1_array_FD)])
xticklabels([1,10,100])
plot_options(gca)
if strcmp(study_case,"simple") && flat_noise && estimate_transient && save_figs
    print(gcf,"report_figures/figures/mean_cost_function_FD_simple_flat","-depsc")
end

TD_keys = "TD_" + l1_array_TD;
minimum = Inf;
for i = 1:numel(TD_keys)
    if mean_Jmr_map(TD_keys(i)) < minimum
        min_TD_key = TD_keys(i);
        minimum = mean_Jmr_map(min_TD_key);
    end
end

FD_keys_true = "FD_true_" + l1_array_FD;
minimum = Inf;
for i = 1:numel(FD_keys_true)
    if mean_Jmr_map(FD_keys_true(i)) < minimum
        min_FD_key_true = FD_keys_true(i);
        minimum = mean_Jmr_map(min_FD_key_true);
    end
end

FD_keys_false = "FD_false_" + l1_array_FD;
minimum = Inf;
for i = 1:numel(FD_keys_false)
    if mean_Jmr_map(FD_keys_false(i)) < minimum
        min_FD_key_false = FD_keys_false(i);
        minimum = mean_Jmr_map(min_FD_key_false);
    end
end

max_l1 = floor(N/2);

summary_keys = [min_TD_key,"TD_" + max_l1,...
    min_FD_key_false,min_FD_key_true,...
    "FD_false_no_l1","FD_true_no_l1",...
    "WNLS_false_LS","WNLS_true_LS",...
    "WNLS_false_actual","WNLS_true_actual"];
if strcmp(study_case,"undermodeled")
    summary_keys = [summary_keys,"actual"];
end
J_summary = getValues(mean_J_map,summary_keys);
Jmr_summary = getValues(mean_Jmr_map,summary_keys);
RMS_summary = getValues(RMS_map,summary_keys);
titles_summary = getTitles(titles_map,summary_keys);

figure
hold on
plot(0.5*db(J_summary),'b--.','Linewidth',1,'Markersize',15,'Displayname','J')
plot(0.5*db(Jmr_summary),'r:o','Linewidth',1,'Markersize',15,'Displayname','J_{mr}')
xticks(1:numel(Jmr_summary))
xticklabels(titles_summary)
legend("Location","Best")
title("Mean cost function at excited frequencies | " + NR_trials + " trials")
set(gca,'fontsize',12)
set(gca,'linewidth',2)
ylabel("dB")
xlim([0.5,numel(Jmr_summary)+0.5])
grid on
set(gca, 'XTickLabelRotation', 45)

indices = [1,4,5];
styles = ["b--","r:","m-.","c-"];
figure
hold on
plot(f_exc/Fs,db(Mk),'k-','Linewidth',2,'Markersize',15,'Displayname','M')
for i = 1:numel(indices)
    plot(f_exc/Fs,db(RMS_summary(:,i)),styles(i),'Linewidth',2,'Markersize',15,...
        'Displayname',"RMS (" + titles_summary(indices(i)) + ")")
end
legend
grid on
ylabel("dB")
xlabel("Normalized frequency")
title("RMS[M - CL] over " + NR_trials + " trials")

estimate_transient = true;
tit1 = "Resulting closed loop system for one noise realization | FD";
tit2 = "Transient";
if ~estimate_transient
    tit2 = tit2 + " not";
end
tit2 = tit2 + " suppressed | ";
if flat_noise
    tit2 = tit2 + "white noise";
else
    tit2 = tit2 + "coloured noise";
end
if estimate_transient
    temp_keys = ["FD_true_1",min_FD_key_true,"FD_true_no_l1"];
else
    temp_keys = ["FD_false_1",min_FD_key_false,"FD_false_no_l1"];
end
names = "l_1 = " + getValues(l1_map,temp_keys);
figure('Units','Normalized','Position',[0.1728    0.4688    0.4876    0.4919])
hold on
plot(f_exc/Fs,db(Mk),'k-','Linewidth',3,'Markersize',15,'Displayname','M')
for i = 1:numel(temp_keys)
    CL = getValues(CL_map,temp_keys(i));
    plot(f_exc/Fs,db(CL(:,1)),styles(i),'Linewidth',3,'Markersize',15,...
        'Displayname',names(i))
end
legend
plot_options(gca)
ylabel("dB")
xlabel("Normalized frequency")
title([tit1,tit2])
if strcmp(study_case,"simple") && flat_noise && estimate_transient && save_figs
    print(gcf,"report_figures/figures/CL_FD_simple_flat","-depsc")
end

%% tables
if strcmp(study_case,"simple") && flat_noise
    table1 = [getValues(mean_J_map,["FD_true_no_l1","FD_false_no_l1"]).',...
            getValues(mean_Jmr_map,["FD_true_no_l1","FD_false_no_l1"]).'];
    table1 = round(0.5*db(table1),2);
    
    keys3 = [min_TD_key,min_FD_key_false,min_FD_key_true,"WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table3 = [getValues(mean_J_map,keys3).',...
            getValues(mean_Jmr_map,keys3).'];
    table3 = round(0.5*db(table3),2);
end

if strcmp(study_case,"simple") && ~flat_noise
    table2 = [getValues(mean_J_map,["FD_true_no_l1","FD_false_no_l1"]).',...
            getValues(mean_Jmr_map,["FD_true_no_l1","FD_false_no_l1"]).'];
    table2 = round(0.5*db(table2),2);
    
    keys4 = [min_TD_key,min_FD_key_false,min_FD_key_true,"WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table4 = [getValues(mean_J_map,keys4).',...
            getValues(mean_Jmr_map,keys4).'];
    table4 = round(0.5*db(table4),2);
end

if strcmp(study_case,"long_transient") && flat_noise    
    keys5 = [min_TD_key,"TD_127",min_FD_key_false,"FD_false_no_l1",min_FD_key_true,"FD_true_no_l1","WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table5 = [getValues(mean_J_map,keys5).',...
            getValues(mean_Jmr_map,keys5).'];
    table5 = round(0.5*db(table5),2);
    
    failed_WNLS = "WNLS_false_LS";
    J = getValues(J_map,failed_WNLS);
    Jmr = getValues(Jmr_map,failed_WNLS);
    failed = J > 1;
    failedmr = Jmr > 1;
    mean_J = mean(J(~failed));
    mean_Jmr = mean(Jmr(~failedmr));
    mean_Jdb = round(0.5*db(mean_J),2);
    mean_Jmrdb = round(0.5*db(mean_Jmr),2);
end

if strcmp(study_case,"long_transient") && ~flat_noise    
    keys6 = [min_TD_key,"TD_127",min_FD_key_false,"FD_false_no_l1",min_FD_key_true,"FD_true_no_l1","WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table6 = [getValues(mean_J_map,keys6).',...
            getValues(mean_Jmr_map,keys6).'];
    table6 = round(0.5*db(table6),2);
    
    failed_WNLS = "WNLS_false_LS";
    J = getValues(J_map,failed_WNLS);
    Jmr = getValues(Jmr_map,failed_WNLS);
    failed = J > 1;
    failedmr = Jmr > 1;
end

if strcmp(study_case,"undermodeled") && flat_noise    
    keys7 = ["actual",min_TD_key,"TD_127",min_FD_key_false,"FD_false_no_l1",min_FD_key_true,"FD_true_no_l1","WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table7 = [getValues(mean_J_map,keys7).',...
            getValues(mean_Jmr_map,keys7).'];
    table7 = round(0.5*db(table7),2);
    
    failed_WNLS = "WNLS_true_LS";
    J = getValues(J_map,failed_WNLS);
    Jmr = getValues(Jmr_map,failed_WNLS);
    failed = J > 0.01 | isnan(J);
    failedmr = Jmr > 0.01 | isnan(Jmr);
    mean_J = mean(J(~failed));
    mean_Jmr = mean(Jmr(~failedmr));
    mean_Jdb = round(0.5*db(mean_J),2);
    mean_Jmrdb = round(0.5*db(mean_Jmr),2);
    
    failed_WNLS2 = "WNLS_true_actual";
    J2 = getValues(J_map,failed_WNLS2);
    Jmr2 = getValues(Jmr_map,failed_WNLS2);
    failed2 = J2 > 0.01;
    failedmr2 = Jmr2 > 0.01;
    mean_J2 = mean(J2(~failed2));
    mean_Jmr2 = mean(Jmr2(~failedmr2));
    mean_Jdb2 = round(0.5*db(mean_J2),2);
    mean_Jmrdb2 = round(0.5*db(mean_Jmr2),2);
end

%%
if strcmp(study_case,"undermodeled") && ~flat_noise    
    keys8 = ["actual",min_TD_key,"TD_127",min_FD_key_false,"FD_false_no_l1",min_FD_key_true,"FD_true_no_l1","WNLS_false_LS","WNLS_false_actual","WNLS_true_LS","WNLS_true_actual"];
    table8 = [getValues(mean_J_map,keys8).',...
            getValues(mean_Jmr_map,keys8).'];
    table8 = round(0.5*db(table8),4);
    
    failed_WNLS = "WNLS_true_LS";
    J = getValues(J_map,failed_WNLS);
    Jmr = getValues(Jmr_map,failed_WNLS);
    failed = J > 0.01 | isnan(J);
    failedmr = Jmr > 0.01 | isnan(Jmr);
    mean_J = mean(J(~failed));
    mean_Jmr = mean(Jmr(~failedmr));
    mean_Jdb = round(0.5*db(mean_J),2);
    mean_Jmrdb = round(0.5*db(mean_Jmr),2);
    
    failed_WNLS2 = "WNLS_true_actual";
    J2 = getValues(J_map,failed_WNLS2);
    Jmr2 = getValues(Jmr_map,failed_WNLS2);
    failed2 = J2 > 0.01 | isnan(J2);
    failedmr2 = Jmr2 > 0.01 | isnan(Jmr2);
    mean_J2 = mean(J2(~failed2));
    mean_Jmr2 = mean(Jmr2(~failedmr2));
    mean_Jdb2 = round(0.5*db(mean_J2),2);
    mean_Jmrdb2 = round(0.5*db(mean_Jmr2),2);
end

%% functions
function vals = getValues(map,keys)
    vals = cell2mat(values(map,cellstr(keys)));
end

function vals = getTitles(titles_map,keys)
    vals = string(values(titles_map,cellstr(keys)));
end

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end