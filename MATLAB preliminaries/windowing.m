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

Gest = Y(k_exc+1,:)./U(k_exc+1);

w = sin(pi*n/N).^2;
W = 1/N*fft(w);

uw = u.*w;
yw = y.*w;

Uw = fft(uw);
Yw = fft(yw);
Gest_windowed = Yw(k_exc+1,:)./Uw(k_exc+1);


fexc = f(k_exc+1);
G_actual_exc = freqresp(Gz,fexc,'Hz');
G_actual_exc = squeeze(G_actual_exc);


figure
hold on
plot(fall/Fs,db(Gall),'Linewidth',2,'Displayname',"|G|")
plot(fexc/Fs,db(G_actual_exc-Gest),'r.--','Markersize',15,'Displayname',"|G_{est}-G|")
plot(fexc/Fs,db(G_actual_exc-Gest_windowed),'ms:','Markersize',12,'Displayname',"|G_{est}-G| (Hann window)")
xlabel('Normalized frequency')
plot_options(gca)
title("Magnitude (dB)")
legend('Location','Best')
print(gcf,'figures/hann_window','-depsc')

function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end