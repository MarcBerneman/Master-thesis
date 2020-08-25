clear
close all

load('non-parametric_est_P=2_M=20_RMS=100mV')
Fs = 78.125e3;
N = 2048;
f = (0:N-1)*Fs/N;
fmax = 10e3;
kmax = floor(fmax/Fs*N);
ExcitedHarm = 1:kmax;

f = f(ExcitedHarm+1);
GBLA = GBLA(ExcitedHarm+1);

%% Find controller
% Controller
s = tf('s');
beta(1,1) = tf(1);
beta(2,1) = 1/s;
beta(3,1) = s;
beta(4,1) = s^2;
% beta(5,1) = 1/s^2;


w01 = 2*pi*5000;
zeta1 = 1;
M1 = w01^2/(s^2 + 2*zeta1*w01*s + w01^2);
w02 = 2*pi*8000;
zeta2 = 1;
M2 = w02^2/(s^2 + 2*zeta2*w02*s + w02^2);
M = M1*M2;

save("M","M")

figure
step(M,(0:100)*0.00001)

k_all = (0:N-1).';
f_all = k_all*Fs/N;
M_all = squeeze(freqresp(M,f_all,'Hz'));
beta_all = squeeze(freqresp(beta,f_all,'Hz')).';
F_all = ones(N,1);
F_all(abs(f_all) > 4000) = 0;
FD_method = FD_controller(M_all,F_all,beta_all);
[rho,Delta] = FD_method.fast_optimize_no_l1(GBLA,ExcitedHarm);
% [rho,Delta] = FD_method.optimize_no_l1(GBLA,ExcitedHarm,1);
K_all = beta_all*rho;
K = beta.'*rho;

save("rho","rho")
save("K","K")

Mk = M_all(ExcitedHarm+1);
Kk = K_all(ExcitedHarm+1);
CLk = Kk.*GBLA./(1+Kk.*GBLA);
Fk = F_all(ExcitedHarm+1);
save("CLk","CLk")
save("Kk","Kk")

%% Figures
figure
hold on
plot(f/1000,db(GBLA),'b-','linewidth',3,'Displayname',"$\mathbf{G_{\mathbf{\mathrm{BLA}}}}$")
plot(f/1000,db(Mk),'r--','linewidth',3,'Displayname',"$\mathbf{M}$")
plot(f/1000,db(CLk),'g-.','linewidth',3,'Displayname',"$\mathbf{CL}$")
plot(f/1000,db(Delta),'k:','linewidth',3,'Displayname',"$\mathbf{M - (1-M)K(\rho)G}$")
xlabel('f (kHz)')
ylabel('Magnitude (dB)')
plot_options(gca)
legend('Location','Best','Interpreter','Latex','fontsize',14)
print("real_system_nonparametric_closed_loop",'-depsc')

figure
hold on
% plot(f/1000,db(GBLA),'b-','linewidth',3,'Displayname',"$\mathbf{G_{\mathbf{\mathrm{BLA}}}}$")
plot(f/1000,db(Mk),'r--','linewidth',3,'Displayname',"$\mathbf{M}$")
plot(f/1000,db(CLk),'g-.','linewidth',3,'Displayname',"$\mathbf{CL}$")
xlabel('f (kHz)')
ylabel('Magnitude (dB)')
plot_options(gca)
legend('Location','Best','Interpreter','Latex','fontsize',14)
print("presentation_real_system_nonparametric_closed_loop",'-dpng','-r300')

%%
figure
hold on
plot(f/1000,db(GBLA),'b-','linewidth',2,'Displayname',"G_{BLA}")
plot(f/1000,db(Fk),'r--','linewidth',2,'Displayname',"F")
xlabel('f (kHz)')
ylabel('Magnitude (dB)')
plot_options(gca)
legend('Location','Best')

%%
[Bm,Am] = tfdata(M,'v');
stability_M = tf(Bm,Am-Bm);
poles = pzmap(stability_M);

figure
plot(poles,'x','Linewidth',4,'Markersize',15)
plot_options(gca)
grid off
xlabel('Real axis (seconds^{-1})')
ylabel('Imaginary axis (seconds^{-1})')
grid on
title("Location of the poles of M/(1-M)")
print("M_one_minus_M_poles",'-depsc')

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end