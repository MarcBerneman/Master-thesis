clear
close all

% file = "data/K_robust_input_P=2_M=20_RMS=100mV";
% file = "data/K_robust_input_P=10_M=100_RMS=100mV";
% file = "data/CL_robust_input_P=20_M=20_RMS=100mV_at_G_in_CL";

file = "data/K_robust_input_P=2_M=20_RMS=100mV_at_K";
% file = "data/CL_robust_input_P=2_M=20_RMS=100mV_at_G_in_CL";
% file = "data/G_robust_input_P=2_M=20_RMS=100mV";

load(file)
% file = "data/K_robust_input_P=2_M=20_RMS=100mV_at_K";
N = Npp;
M = size(u,2);
RMS = round(rms(u(:,1)),4);

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
% method.dof = 60;
% method.order = 2; % order of transient
[CZ, Z, freq, GBLA, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
GBLA = squeeze(GBLA);
sigmaBLA = squeeze(CvecG.NL);
sigmaBLAn = squeeze(CvecG.n);
sigmaGs = squeeze(CvecG.S);

% R = 2;
% q = 1;
% n_rec = ceil((q + R+1)/2);
% dof_rec = 2*n_rec - (R + 1);

%% figures
load('CLk')
load('Kk')
figure
hold on

if strcmp(file,"data/K_robust_input_P=2_M=20_RMS=100mV_at_K")
    plot(fexc,db(Kk),'m:','Linewidth',2,'Displayname','|K|^2 (optimal)')
    plot(fexc,db(GBLA),'k-','Linewidth',2,'Displayname','|K|^2 (measured)') %20 log(|G|) = 10 log(|G|^2)
elseif strcmp(file,"data/CL_robust_input_P=2_M=20_RMS=100mV_at_G_in_CL")
    plot(fexc,db(CLk),'m:','Linewidth',2,'Displayname','|CL|^2 (optimal)')
    plot(fexc,db(GBLA),'k-','Linewidth',2,'Displayname','|CL|^2 (measured)') %20 log(|G|) = 10 log(|G|^2)
else
    plot(fexc,db(GBLA),'k-','Linewidth',2,'Displayname','|G_{BLA}|^2') %20 log(|G|) = 10 log(|G|^2)
    plot(fexc,0.5*db(sigmaBLA),'b.','Displayname','Total variance') %10 log(variance)
    plot(fexc,0.5*db(sigmaBLAn),'ro','Displayname','Noise variance')
    plot(fexc,0.5*db(sigmaGs),'gx','Displayname','NL distortion 1 realization')
end
grid on
title({"Robust method",...
     "N = " + N + " | P = " + P + " | M = " + M + " | RMS = " + RMS*1000 + "mV | f_{max} = " + fmax/1000 + "kHz"})
legend('Location','Best','fontsize',10)
xlim([0 fmax])
ylabel("dB")
plot_options(gca)
if strcmp(file,"data/G_robust_input_P=2_M=20_RMS=100mV")
    xlabel('f (kHz)')
    xticklabels(xticks/1000)
    print(gcf,'robust_method_G.eps','-depsc')
elseif strcmp(file,"data/K_robust_input_P=2_M=20_RMS=100mV_at_K")
    xlabel('f (Hz)')
    set(gca,'xscale','log')
    print(gcf,'robust_method_K.eps','-depsc')
elseif strcmp(file,"data/CL_robust_input_P=2_M=20_RMS=100mV_at_G_in_CL")
    xlabel('f (Hz)')
    print(gcf,'robust_method_CL.eps','-depsc')
end


%%
figure
for i = 2:-1:1
    subplot(2,1,3-i)
    hold on
    plot(fexc,db(squeeze(Z(i,1,:))),'k-','Linewidth',1.5,'Displayname','Mean')
    plot(fexc,0.5*db(squeeze(CZ.NL(i,i,1,:))),'b.','Displayname','Total variance')
    plot(fexc,0.5*db(squeeze(CZ.n(i,i,1,:))),'ro','Displayname','Noise variance')
    plot(fexc,0.5*db(squeeze(CZ.S(i,i,1,:))),'gx','Displayname','NL distortion 1 realization')
    hold off
    xticklabels(xticks/1000)
    xlabel('f (kHz)')
    grid on
    legend
    ylabel("dB")
    set(gca,'fontsize',10)
    set(gca,'linewidth',1.5)
%     set(gca,'xscale','log')
%     xlim([min(fexc),max(fexc)])
end
subplot(2,1,1)
title("Input")
subplot(2,1,2)
title("Output")


function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end