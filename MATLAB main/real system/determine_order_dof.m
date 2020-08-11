clear
close all

% file = "data/K_robust_input_P=2_M=20_RMS=100mV";
% file = "data/K_robust_input_P=10_M=100_RMS=100mV";
% file = "data/K_robust_input_P=2_M=20_RMS=100mV_at_K";
% file = "data/CL_robust_input_P=20_M=20_RMS=100mV_at_G_in_CL";
file = "data/G_robust_input_P=2_M=20_RMS=100mV";
load(file)

N = Npp;
M = size(u,2);
RMS = 0.1;

fres = fs/N;
f = (0:N-1)*fres;
fmax = 10e3;
kmax = floor(fmax/fres);
kexc = 1:kmax;
fexc = kexc*fres;
F = numel(kexc);

%% robust method
ny = 1; nu = 1;
data = struct('y', zeros(ny,nu,M,P*N), 'u', zeros(nu,nu,M,P*N) , 'r', zeros(nu,nu,M,N));
data.y(1,1,:,:) = y_meas(:,:).';
data.u(1,1,:,:) = u_meas(:,:).';
data.r(1,1,:,:) = u(:,:).';
data.Ts = 1/fs;
data.N = N;
data.ExcitedHarm = kexc;
method.transient = 1; % take transient into account

orders = 0:2:8;
for i = numel(orders):-1:1
    method.order = orders(i); % order of transient
    [CZ, Z, freq, GBLA, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
    sigmaBLA_orders{i} = squeeze(CvecG.NL);
end

method.order = 2;
dofs = 0:10;
for i = numel(dofs):-1:1
    method.dof = dofs(i); % order of transient
    [CZ, Z, freq, GBLA, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
    sigmaBLA_dofs{i} = squeeze(CvecG.NL);
end
%% figures
figure
hold on
for i = 1:numel(orders)
    plot(fexc,0.5*db(sigmaBLA_orders{i}),'Displayname',"R = " + orders(i)) %10 log(variance)
end
xticklabels(xticks/1000)
xlabel('f (kHz)')
grid on
title({"Robust method | Total variance",...
     "N = " + N + " | P = " + P + " | M = " + M + " | RMS = " + RMS*1000 + "mV | f_{max} = " + fmax/1000 + "kHz"})
legend('fontsize',10)
xlim([0 fmax])
ylabel("dB")
set(gca,'fontsize',10)
set(gca,'linewidth',1.5)

figure
hold on
for i = 1:numel(dofs)
    plot(fexc,0.5*db(sigmaBLA_dofs{i}),'Displayname',"dof = " + dofs(i)) %10 log(variance)
end
xticklabels(xticks/1000)
xlabel('f (kHz)')
grid on
title({"Robust method | Total variance",...
     "N = " + N + " | P = " + P + " | M = " + M + " | RMS = " + RMS*1000 + "mV | f_{max} = " + fmax/1000 + "kHz"})
legend('fontsize',10)
xlim([0 fmax])
ylabel("dB")
set(gca,'fontsize',10)
set(gca,'linewidth',1.5)