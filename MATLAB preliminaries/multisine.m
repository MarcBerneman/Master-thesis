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
n = 0:N-1;
t = (0:N-1)*Ts;
f = (0:N-1)*Fs/N;

P = 20;
u_rep = repmat(u,P,1);
y_rep = lsim(Gz,u_rep,(0:N*P-1)*Ts);
y = y_rep((1:N)+(P-1)*N);
Y = 1/N*fft(y);

Gest = Y(k_exc+1)./U(k_exc+1);
fexc = f(k_exc+1);

figure
subplot(2,1,1)
plot(n,u,'-','Linewidth',2)
xlabel('n')
title("u(n)")
plot_options(gca)
subplot(2,1,2)
stem(f(1:N/2)/Fs,abs(U(1:N/2)),'Linewidth',2)
xlabel('Normalized frequency')
title("|U|")
plot_options(gca)
print(gcf,'figures/MS_u','-depsc')

subplot(2,1,1)
xlabel('Time index')
ylabel('Arbitrary unit')
subplot(2,1,2)
ylabel('Arbitrary unit')
print(gcf,'figures/presentation_MS_u','-dpng','-r300')

figure
subplot(2,1,1)
plot(n,y,'-','Linewidth',2)
xlabel('n')
title("y(n)")
plot_options(gca)
subplot(2,1,2)
plot(f(1:N/2)/Fs,db(Y(1:N/2)),'.','Color',[0 0.4470 0.7410],'Markersize',15,'Linewidth',2)
xlabel('Normalized frequency')
title("|Y|_{dB}")
plot_options(gca)
ylim([-40,20])
xlim([min(f(1:N/2)),max(f(1:N/2))]/Fs)
print(gcf,'figures/MS_y','-depsc')

subplot(2,1,1)
xlabel('Time index')
ylabel('Arbitrary unit')
subplot(2,1,2)
ylabel('Arbitrary unit')
print(gcf,'figures/presentation_MS_y','-dpng','-r300')


figure
subplot(2,1,1)
plot(f(1:N/2)/Fs,db(Y(1:N/2)),'.','Color',[0 0.4470 0.7410],'Markersize',15,'Linewidth',2)
title("|Y|_{dB}")
ylabel('Arbitrary unit')
plot_options(gca)
ylim([-40,20])
xlim([min(f(1:N/2)),max(f(1:N/2))]/Fs)
subplot(2,1,2)
stem(f(1:N/2)/Fs,abs(U(1:N/2)),'Linewidth',2)
xlabel('Normalized frequency')
title("|U|")
ylabel('Arbitrary unit')
plot_options(gca)
print(gcf,'figures/presentation_MS_YU','-dpng','-r300')

figure
hold on
plot(fall/Fs,db(Gall),'Linewidth',2)
plot(fexc/Fs,db(Gest),'r.','Markersize',15)
xlabel('Normalized frequency')
plot_options(gca)
title("|G|_{dB}")
legend("Actual","Estimated")
print(gcf,'figures/MS_G','-depsc')

grid off
legend off
xticks([])
yticks([])
xlabel("")
set(gca,'XColor','w','YColor','w')
title("")
print(gcf,'figures/presentation_freq_dom_title','-dpng','-r300')

function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end