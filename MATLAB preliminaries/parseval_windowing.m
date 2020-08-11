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

energy_FD0 = H0'*H0*N;
energy_FD = mean(sum(abs(H).^2)*N,2);
energy_TD0 = cumsum(h0.^2);
energy_TD = mean(cumsum(h.^2),2);

figure
hold on
stairs(n,h0,'b-','Linewidth',2,'Displayname','Noiseless')
stairs(n,h,'r:','Linewidth',2,'Displayname',"Noisy (\sigma_h = " + sigma_h + ")")
legend('Location','Best')
xlabel("Samples n")
title("Impulse response of H(z^{-1})")
grid on
plot_options(gca)

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

l1 = 10;
w = zeros(N,1);
w(1:l1+1) = 1;
W = 1/N*fft(w);

h_cutoff = h.*w;
h_cutoff(l1+2:end) = 0;
H_cutoff = 1/N*fft(h_cutoff);
H_cutoff2 = cconv(H,W,N);

figure
subplot(2,2,1)
stairs(n,h)
subplot(2,2,2)
stem(n,abs(H))
subplot(2,2,3)
stairs(n,h_cutoff)
subplot(2,2,4)
hold on
stem(n,abs(H_cutoff),'o')
stem(n,abs(H_cutoff2),'x')
%% functions
function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end
