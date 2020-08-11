clear
close all

w0 = 2*pi*0.1;
zeta = 0.1;

s = tf('s');
G = 0.5*w0^2/(s^2 + 2*zeta*w0*s + w0^2);
Ts = 1;
Fs = 1/Ts;

% Gz = c2d(G,Ts,'zoh');
% Gz = tf([1,-3.22,4.59,-3.22,1],[1,-0.94,0.08],Ts,'Variable','z^-1');
Gz = zpk([exp(2*pi*0.1j),exp(-2*pi*0.1j),exp(2*pi*0.25j),exp(-2*pi*0.25j),-1],[0,0,0.35*exp(-2*pi*0.4j),0.35*exp(2*pi*0.4j),0,0.9],1,[]);
Gz = Gz/dcgain(Gz);

figure
pzmap(Gz)

Su = Gz;
Sy = Gz;
save("Noise_filters.mat",'Su','Sy')

N = 255;
if even(N)
    k_max = N/2-1;
else
    k_max = (N-1)/2;
end
f = (0:k_max)*Fs/N;

figure
impulse(Su)
title("Impulse response of S_u")

Su = squeeze(freqresp(Su,f,'Hz'));
Sy = squeeze(freqresp(Sy,f,'Hz'));

figure
hold on
plot(f,db(Su),'b-','Linewidth',2,'Displayname',"|S_u|")
plot(f,db(Sy),'r:','Linewidth',2,'Displayname',"|S_y|")
xlabel('Frequency (Hz)')
plot_options(gca)
title("Magnitude (dB)")
legend('Location','Best')


function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end