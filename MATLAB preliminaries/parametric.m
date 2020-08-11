clear
close all
rng(0)

w0 = 2*pi*0.3;
zeta = 0.01;

s = tf('s');
G = w0^2/(s^2 + 2*zeta*w0*s + w0^2);
Fs = 2;
Ts = 1/Fs;

Gz = c2d(G,Ts,'zoh');

fall = linspace(0,Fs/2,10000);
Gall = freqresp(Gz,fall,'Hz');
Gall = squeeze(Gall);

N = 100;
fmin = 0.1;
fmax = 0.7;
df = 0.02;
RMS = 1;
[u,k_exc] = generate_MS(N,Fs,RMS,fmin,fmax,df);
U = 1/N*fft(u);
n = (0:N-1)';
t = n*Ts;
f = n*Fs/N;

y = lsim(Gz,u,t);
Y = 1/N*fft(y);

Y = Y(k_exc+1);
U = U(k_exc+1);

na = 2;
nb = 2;
nI = max(na,nb)-1;

Omega_k = exp(-1j*2*pi*f(k_exc+1)*Ts);
Ka = Y.*Omega_k.^(0:na);
Kb = -U.*Omega_k.^(0:nb);
KI = -Omega_k.^(0:nI);
K = [Ka,Kb,KI];
rho = null(K);
rho = rho/rho(1);

a = real(rho(1:na+1)).';
b = real(rho(na+2:na+nb+2)).';
Gest = tf(b,a,Ts,'Variable','z^-1');
Gest_all = squeeze(freqresp(Gest,fall,'Hz'));

% K2 = [Ka,Kb];
% rho2 = null(K2);
% rho2 = rho2/rho2(1);
% 
% a2 = real(rho2(1:na+1)).';
% b2 = real(rho2(na+2:na+nb+2)).';
% Gest2 = tf(b2,a2,Ts,'Variable','z^-1');
% Gest_all2 = squeeze(freqresp(Gest2,fall,'Hz'));

figure
hold on
plot(fall/Fs,db(Gall),'Linewidth',2,'Displayname',"|G|")
plot(fall/Fs,db(Gall-Gest_all),'r--','Linewidth',2,'Displayname',"|G_{est}-G|")
% plot(fexc/Fs,db(G_actual_exc-Gest_windowed),'ms:','Markersize',12,'Displayname',"|G_{est}-G| (Hann window)")
xlabel('Normalized frequency')
plot_options(gca)
title("Magnitude (dB)")
legend('Location','Best')
print(gcf,'figures/parametric_transient','-depsc')

function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end