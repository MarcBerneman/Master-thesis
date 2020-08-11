clear
close all
save_fig = true;

rng_index = 0;
rng(rng_index)

%% load configurations
study_case = "simple";
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
estimate_transient = true;
NR_trials = 100;
% sigma_y = 0.05;

load(study_case + "/System_data")
Fs = 1/Ts;

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
d = randn(Ntot,1)*sigma_d;
u0 = r0 + d;
y0 = lsim(G,u0,t);

% Noise filters
if ~flat_noise
    load("Noise_filters")
end

% Noise
ey = randn(numel(u0),NR_trials);
eu = randn(numel(u0),NR_trials);
ey = ey.*sigma_y;
eu = eu.*sigma_u;

if flat_noise
    vy = ey;
    vu = eu;
else
    for i = NR_trials:-1:1
        vy(:,i) = lsim(Sy,ey(:,i),0:Ntot-1);
        vu(:,i) = lsim(Su,eu(:,i),0:Ntot-1);
    end
end
y = y0 + vy;
u = u0 + vu;

% Transient removal
r0 = r0((P_SS+P_zero_input)*N+1:end);
u0 = u0((P_SS+P_zero_input)*N+1:end);
y0 = y0((P_SS+P_zero_input)*N+1:end);
vy = vy((P_SS+P_zero_input)*N+1:end,:);
y = y((P_SS+P_zero_input)*N+1:end,:);
vu = vu((P_SS+P_zero_input)*N+1:end,:);
u = u((P_SS+P_zero_input)*N+1:end,:);
t = t(1:end-(P_SS+P_zero_input)*N);
Ntot = Ntot - (P_SS+P_zero_input)*N;

k_all = (0:N-1).';
f_all = k_all*Fs/N;

% nonparametric estimate of the open loop plant
order_array = 2:2:8;
for i = numel(order_array):-1:1
    order = order_array(i);
    dof = 1;
    [G_est{i},varG{i}] = nonparametric_estimate(r0,u(:,1),y(:,1),P,Ts,ExcitedHarm,order,dof,estimate_transient);
end
figure
hold on
for i = 1:numel(order_array)
    plot(f_plot/Fs,0.5*db(varG{i}),'Linewidth',2)
end
legend("R = " + order_array)
title("\sigma_G")

grid on

G_actual = squeeze(freqresp(G,f_plot,'Hz'));
figure
hold on
for i = 1:numel(order_array)
    plot(f_plot/Fs,db(G_est{i}-G_actual),'Linewidth',2)
end
title("|G_{est}-G_{actual}|^2_{dB}")
legend("R = " + order_array)
grid on

order = 6;
dof_array = 0:5;
for i = numel(dof_array):-1:1
    dof = dof_array(i);
    [G_est{i},varG{i}] = nonparametric_estimate(u0,u(:,1),y(:,1),P,Ts,ExcitedHarm,order,dof,estimate_transient);
end
figure
hold on
for i = 1:numel(dof_array)
    plot(f_plot/Fs,0.5*db(varG{i}),'Linewidth',2)
end
title("\sigma_G")
legend("dof = " + dof_array)
grid on

G_actual = squeeze(freqresp(G,f_plot,'Hz'));
figure
hold on
for i = 1:numel(dof_array)
    plot(f_plot/Fs,db(G_est{i}-G_actual),'Linewidth',2)
end
title("|G_{est}-G_{actual}|^2_{dB}")
legend("dof = " + dof_array)
grid on

order = 6;
dof = 6;
RMS_est = zeros(size(G_actual));
RMS_est_no_transient = zeros(size(G_actual));
for i = NR_trials:-1:1
    [G_est,varG] = nonparametric_estimate(r0,u(:,i),y(:,i),P,Ts,ExcitedHarm,order,dof,estimate_transient);
    [G_est_no_transient,varG_no_transient] = nonparametric_estimate(r0,u(:,i),y(:,i),P,Ts,ExcitedHarm,order,dof,false);
    RMS_est = RMS_est + abs(G_est - G_actual).^2;
    RMS_est_no_transient = RMS_est_no_transient + abs(G_est_no_transient - G_actual).^2;
end
RMS_est = sqrt(RMS_est/NR_trials);
RMS_est_no_transient = sqrt(RMS_est_no_transient/NR_trials);

figure
hold on
plot(f_plot/Fs,db(G_est),'k--.','Markersize',15,'Linewidth',2,'Displayname',"G_{est}")
plot(f_plot/Fs,0.5*db(varG),'b--','Linewidth',2,'Displayname',"\sigma_G")
grid on
legend


figure('Units','Normalized','Position',[0.1728    0.2688    0.4876    0.4919])
hold on
plot(f_plot/Fs,db(RMS_est),'b-','Linewidth',3)
plot(f_plot/Fs,db(RMS_est_no_transient),'r--','Linewidth',3)
xlabel("Normalized frequency")
tit1 = "MSE[G_{est}] (dB) (" + NR_trials + " realizations)";
legend(["Transient removed","No transient removed"])
plot_options(gca)
if strcmp(study_case,"simple") && ~flat_noise && save_fig
    tit2 = "Simple system | coloured noise";
    title([tit1,tit2])
    print("report_figures/figures/MSE_Gest_simple_coloured","-depsc")
else
    title(tit1)
end


figure
hold on
plot(f_plot/Fs,0.5*db(varG),'b-','Linewidth',2)
plot(f_plot/Fs,0.5*db(varG_no_transient),'r--','Linewidth',2)
title("\sigma_G (dB)")
legend(["Transient removed","No transient removed"])
grid on

%% functions
function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end