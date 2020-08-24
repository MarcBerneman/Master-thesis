clear
close all

%% load data
study_case = "undermodeled";
load("../" + study_case + "/System_data")
if ~exist("M","Var")
    M = M_exact;
end

N = 255;
if even(N)
    k_max = N/2-1;
else
    k_max = (N-1)/2;
end
ExcitedHarm = 1:k_max;
Fs = 1/M.Ts;
f = (0:k_max)*Fs/N;
f_exc = f(ExcitedHarm+1);
Gk = squeeze(freqresp(G,f_exc,'Hz'));
Mk = squeeze(freqresp(M,f_exc,'Hz'));

n_rho = 4;
qinv = tf([0,1],1,Ts,"variable","q^-1");
for i = n_rho:-1:0
    beta(i+1,1) = qinv^i / (1-qinv);
end
beta_all = squeeze(freqresp(beta,f_exc,'Hz'));
if n_rho > 0
    beta_all = beta_all.';
end
load("../" + study_case + "/rho_actual")
Kopt = beta_all*rho_actual;
CLopt = Gk.*Kopt./(1+Gk.*Kopt);

%% figures
figure
hold on
plot(f_exc/Fs,db(CLopt),'b-','Linewidth',4,'Markersize',15,'Displayname','CL')
plot(f_exc/Fs,db(Mk),'r:','Linewidth',4,'Markersize',15,'Displayname','M')
legend("Location","Best")
grid on
xlabel("Normalized frequency")
title("Magnitude (dB)")
plot_options(gca)
print(gcf,"figures/undermodeled_optimal","-depsc")

set(gcf,'Units','Normalized','Position',[0.1,0.1,0.8,0.7])
grid off
title("")
legend off
xlabel("")
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'XColor','w','YColor','w')
print(gcf,"figures/presentation_title_page","-djpeg")


function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',18)
    grid on
end