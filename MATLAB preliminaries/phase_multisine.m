clear
close all
rng(0);

N = 100;
Fs = 1;
Ts = 1/Fs;
fmin = 0.01;
fmax = 0.1;
df = 0.01;
RMS = 1;
[u,k_exc] = generate_MS(N,Fs,RMS,fmin,fmax,df);
U = 1/N*fft(u);
[u0,~] = generate_MS(N,Fs,RMS,fmin,fmax,df,0);
U0 = 1/N*fft(u0);
t = (0:N-1)*Ts;
f = (0:N-1)*Fs/N;


figure
subplot(2,1,1)
hold on
plot(t,u0,'b-','Linewidth',2)
plot(t,u,'r:','Linewidth',2)
legend("Linear phase","Random phase")
xlabel('t(s)')
title("u(t)")
plot_options(gca)
subplot(2,1,2)
hold on
stem(f(1:N/2),abs(U0(1:N/2)),'bo','Linewidth',2)
stem(f(1:N/2),abs(U(1:N/2)),'rx','Linewidth',1)
legend("Linear phase","Random phase")
% xlim([0,2*fmax])
xlabel('Frequency (Hz)')
title("|U|")
plot_options(gca)
print(gcf,'figures/MS_u0_u','-depsc')


function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end