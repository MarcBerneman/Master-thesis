clear
close all
rng(0)

N = 51;
n = 0:N-1;
Ts = 1; Fs = 1/Ts;

z = tf('z',Ts);
Gz = (1-0.8)/(z-0.8);

M = 1;
sigma_h = 0.01;
nh = sigma_h*randn(N,M);
h0 = impulse(Gz,n*Ts);
h = h0 + nh;
H0 = 1/N*fft(h0);
H = 1/N*fft(h);

energy_FD0 = 1/N*(H0'*H0);
energy_FD = 1/N*sum(abs(H).^2);
energy_TD0 = 1/N^2*cumsum(h0.^2);
energy_TD = 1/N^2*cumsum(h.^2);

figure
hold on
stairs(n,h0,'b-','Linewidth',2,'Displayname','Noiseless')
stairs(n,h,'r:','Linewidth',2,'Displayname',"Noisy (\sigma_h = " + sigma_h + ")")
legend('Location','Best')
xlabel("Samples n")
title("Impulse response of H(z^{-1})")
grid on
plot_options(gca)
print(gcf,'figures/parseval_signal','-depsc')

figure
hold on
line2 = yline(energy_FD,':','Color',[1,0,0],'Linewidth',5,'Displayname',"FD noisy (\sigma_h = " + sigma_h + ")");
line2.Alpha = 0.3;
plot(n,energy_TD,'r:','Linewidth',2,'Displayname',"TD noisy (\sigma_h = " + sigma_h + ")")
line = yline(energy_FD0,'-.','Color',[0,0,1],'Linewidth',5,'Displayname',"FD noiseless");
line.Alpha = 0.3;
plot(n,energy_TD0,'b-','Linewidth',2,'Displayname','TD noiseless')
title("Energy in the time domain (TD) and frequency domain (FD)")
grid on
xlabel("l_1")
legend('Location','Best')
plot_options(gca)
yticklabels("0.0000" + yticks*1e5)
print(gcf,'figures/parseval_energy','-depsc')

%% functions
function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end
