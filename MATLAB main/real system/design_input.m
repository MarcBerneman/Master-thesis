clear
close all


load("inputs/robust_input_P=2_M=20_RMS=100mV")
Fs = 78.125e3;
Ts = 1/Fs;
N = 2048;
fmax = 10e3;
kmax = floor(fmax/Fs*N);
ExcitedHarm = 1:kmax;

k = 0:N-1;
t = k*Ts;
f = k*Fs/N;
f = f(ExcitedHarm+1);
U = 1/N*fft(u);

%% input for open loop
m = 1;
um = u(:,m);
Um = U(:,m);
Um = Um(ExcitedHarm+1);
figure
subplot(2,1,1)
plot(t,um)
title('u(t)')
xlabel('t (s)')
subplot(2,1,2)
plot(f,abs(Um))
title('|U(f)|')
xlabel('f (Hz)')

rms_t = rms(um);
rms_f = sqrt(2)*sqrt(sum(abs(Um).^2));

%% input for closed loop
load('CLk')
load("non-parametric_est_P=2_M=20_RMS=100mV")
GBLA = GBLA(ExcitedHarm+1);
R_to_U = CLk./GBLA;

R_to_U_complete = zeros(N,1);
R_to_U_complete(ExcitedHarm+1) = R_to_U;
R_to_U_complete = fft(2*real(ifft(R_to_U_complete)));

REF_CL = U./R_to_U_complete;
REF_CL(abs(U) < 1e-9) = 0;
ref_CL = N*ifft(REF_CL);
rms_ref_CL = rms(ref_CL(:,m));

save_data(ref_CL,P,'robust_input_P=2_M=20_RMS=100mV_at_G_in_CL')

figure
subplot(2,1,1)
plot(t,ref_CL(:,m))
title('u(t) (CL)')
xlabel('t (s)')
subplot(2,1,2)
plot(f,db(REF_CL(ExcitedHarm+1,m)))
title('|U(f)|_{dB} (CL)')
xlabel('f (Hz)')

%% input for controller measuring
load('Kk')

Kk_complete = zeros(N,1);
Kk_complete(ExcitedHarm+1) = Kk;
Kk_complete = fft(2*real(ifft(Kk_complete)));

REF_K = U./Kk_complete;
REF_K(abs(U) < 1e-9) = 0;
ref_K = N*ifft(REF_K);
rms_ref_K = rms(ref_K(:,m));

save_data(ref_K,P,'robust_input_P=2_M=20_RMS=100mV_at_K')

figure
% subplot(2,1,1)
% plot(t,ref_K(:,m))
% title('u(t) (CL)')
% xlabel('t (s)')
% subplot(2,1,2)
plot(f,db(REF_K(ExcitedHarm+1,m)),'b-','Linewidth',4)
title('|U(f)|^2_{dB} (1 period, 1 realization)')
xlabel('f (Hz)')
plot_options(gca)
set(gca,'xscale','log')
xlim([min(f),max(f)])
print(gcf,'controller_input_for_meas','-depsc')

function save_data(u,P,name)
    save(name,'u','P')
end

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end