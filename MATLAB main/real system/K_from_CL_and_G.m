clear
close all

file = "data/G_robust_input_P=2_M=20_RMS=100mV";
load(file)
N = Npp;
M = size(u,2);

fres = fs/N;
f = (0:N-1)*fres;
fmax = 10e3;
kmax = floor(fmax/fres);
kexc = 1:kmax;
fexc = kexc*fres;
F = numel(kexc);

ny = 1; nu = 1;
data = struct('y', zeros(ny,nu,M,P*N), 'u', zeros(nu,nu,M,P*N) , 'r', zeros(nu,nu,M,N));
data.y(1,1,:,:) = y_meas(:,:).';
data.u(1,1,:,:) = u_meas(:,:).';
data.r(1,1,:,:) = u(:,:).';
data.Ts = 1/fs;
data.N = N;
data.ExcitedHarm = kexc;
method.transient = 1; % take transient into account
% method.dof = 60;
% method.order = 2; % order of transient
[CZ, Z, freq, GBLA, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
G = squeeze(GBLA);

file = "data/CL_robust_input_P=2_M=20_RMS=100mV_at_G_in_CL";
load(file)
N = Npp;
M = size(u,2);
RMS = 0.1;
data = struct('y', zeros(ny,nu,M,P*N), 'u', zeros(nu,nu,M,P*N) , 'r', zeros(nu,nu,M,N));
data.y(1,1,:,:) = y_meas(:,:).';
data.u(1,1,:,:) = u_meas(:,:).';
data.r(1,1,:,:) = u(:,:).';
data.Ts = 1/fs;
data.N = N;
data.ExcitedHarm = kexc;
method.transient = 1; % take transient into account
% method.dof = 60;
% method.order = 2; % order of transient
[CZ, Z, freq, GBLA, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
CL = squeeze(GBLA);


K_from_G_and_CL = CL./G./(1-CL);

%% figures
load('Kk')
% load('CLk')

figure
hold on
plot(fexc,db(Kk),'m:','Linewidth',4,'Displayname','|K|^2 (optimal)')
plot(fexc,db(K_from_G_and_CL),'k-','Linewidth',4,'Displayname','|K|^2 (from G and CL)')
% plot(fexc,db(CLk),'m:','Linewidth',2,'Displayname','|CL|^2 (optimal)')
% plot(fexc,db(CL_fromKG),'k-','Linewidth',2,'Displayname','|CL|^2 (from K, G)') %20 log(|G|) = 10 log(|G|^2)
grid on
title("K estimated indirectly from G and CL")
legend('Location','Best','fontsize',14)
xlim([0 fmax])
ylabel("dB")
plot_options(gca)
print(gcf,'K_from_G_CL.eps','-depsc')


function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end