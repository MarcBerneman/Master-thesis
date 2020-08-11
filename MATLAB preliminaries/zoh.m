clear
close all

w0 = 2*pi*0.3;
zeta = 0.01;

s = tf('s');
G = w0^2/(s^2 + 2*zeta*w0*s + w0^2);
Fs = 2;
Ts = 1/Fs;

Gz = c2d(G,Ts,'zoh');

fall = linspace(0,Fs,10000);
Gall = squeeze(freqresp(G,fall,'Hz'));
Gzall = squeeze(freqresp(Gz,fall,'Hz'));

figure
hold on
plot(fall,db(Gall),'b-','Linewidth',2,'Displayname',"|G(s)|")
plot(fall,db(Gzall),'r:','Linewidth',2,'Displayname',"|G_{ZOH}(z^{-1})|")
xlabel('Frequency (Hz)')
plot_options(gca)
title("Magnitude (dB)")
legend('Location','Best')
print(gcf,'figures/ZOH','-depsc')

function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end