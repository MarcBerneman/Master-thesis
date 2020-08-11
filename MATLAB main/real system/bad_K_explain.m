clear
close all


file = "data/K_robust_input_P=2_M=20_RMS=100mV_at_K";
% file = "data/CL_robust_input_P=2_M=20_RMS=100mV_at_G_in_CL";
load(file)

N = Npp;
M = size(u,2);
RMS = round(rms(u(:,1)),4);

fres = fs/N;
f = (0:N-1)*fres;
f_full = (0:N*P-1)*fs/(N*P);
fmax = 10e3;
kmax = floor(fmax/fres);
kexc = 1:kmax;
fexc = kexc*fres;
F = numel(kexc);

U = 1/N*fft(u(:,1));
U_meas = 1/(N*P)*fft(u_meas(:,1));
Y_meas = 1/(N*P)*fft(y_meas(:,1));


figure
subplot(2,1,1)
hold on
plot(f(1:end/2)/1000,abs(U(1:end/2)),'linewidth',2,'Displayname','Reference')
plot(f_full(1:end/2)/1000,abs(U_meas(1:end/2)),'.','Displayname',"Measured input (" + P + " periods)")
% ylim([-150,-20])
legend('location','best')
ylabel('Magnitude (dB)')
title('Measurement of the controller for 1 realization')
plot_options(gca)
subplot(2,1,2)
plot(f_full(1:end/2)/1000,db(Y_meas(1:end/2)),'.','Displayname',"Measured output (" + P + " periods)")
legend('location','best')
xlabel('Frequency (kHz)')
ylabel('Magnitude (dB)')
plot_options(gca)
print(gcf,'bad_K_explain','-djpeg')

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end