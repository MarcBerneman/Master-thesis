clear
close all
rng(0)

order = 2;
dof = 9;
windowsize = (dof + order + 1)/2;
sigma_u = 0;
sigma_y = 0.20;

w0 = 2*pi*0.3;
zeta = 0.01;
s = tf('s');
G = w0^2/(s^2 + 2*zeta*w0*s + w0^2);
Fs = 2;
Ts = 1/Fs;
z = tf('z',Ts);
Gz = c2d(G,Ts,'zoh');
N = 2^14;
fmin = Fs/N;
fmax = 5000*Fs/N;
df = Fs/N;

% Fs = 1e3; Ts = 1/Fs;
% z = tf('z',Ts);
% Gz = minreal((0.1925*z^-1 + 0.1884*z^-2)/(1-1.5574*z^-1+0.9382*z^-2));
% N = 2^14;
% fmin = Fs/N;
% fmax = 4096*Fs/N;
% df = Fs/N;

fall = linspace(0,Fs/2,10000);
Gall = freqresp(Gz,fall,'Hz');
Gall = squeeze(Gall);


RMS = 1;
[u0,k_exc] = generate_MS(N,Fs,RMS,fmin,fmax,df);

P = 2;
n = (0:N*P-1)';
u0 = repmat(u0,P,1);
% [bBJ,aBJ]  = butter(2,2*0.133);
Hz = minreal((1-3.2193*z^-1+4.5893*z^-2-3.2193*z^-3+z^-4)/(1-0.9355*z^-1+0.0780*z^-2));
[bBJ,aBJ] = tfdata(Hz,'v');

w = sin(pi*n/(N*P)).^2;
W = fft(w);
NR_trials = 100;
for trial = NR_trials:-1:1
    eu = sigma_u*randn(size(n));
    ey = sigma_y*randn(size(n));
    nu = filter(bBJ,aBJ,eu);
    ny = filter(bBJ,aBJ,ey);

    t = n*Ts;
    f = n*Fs/N;
    u = u0 + nu;
%     y = lsim(Gz,u,t) + ny; %use filtered white noise
    y = lsim(Gz,u,t) + ey; %use white noise
    U = 1/N*fft(u);
    Y = 1/N*fft(y);

    G_est_normal(:,trial) = Y(P*k_exc+1)./U(P*k_exc+1);
    [G_est_robust(:,trial),varG(:,trial)] = nonparametric_estimate(u,u,y,P,Ts,k_exc,order,dof);

    uw = u.*w;
    yw = y.*w;

    Uw = 1/N*fft(uw);
    Yw = 1/N*fft(yw);
    G_est_windowed(:,trial) = Yw(P*k_exc+1,:)./Uw(P*k_exc+1);
end

fexc = f(k_exc+1);
G_actual_exc = squeeze(freqresp(Gz,fexc,'Hz'));
RMS_normal = rms(abs(G_actual_exc-G_est_normal),2);
RMS_windowed = rms(abs(G_actual_exc-G_est_windowed),2);
RMS_robust = rms(abs(G_actual_exc-G_est_robust),2);

%% SNR calculation
y0 = lsim(Gz,u0,t);
SNR = rms(y0)/sigma_y;
SNRdb = db(SNR);

%% figures
diffRMS = db(RMS_robust)-db(RMS_windowed);
mean_window = 50;
fexc_mean = window_mean(fexc,mean_window);
diffRMS_mean = window_mean(diffRMS,mean_window);
RMS_normal_mean = window_mean(RMS_normal,mean_window);
RMS_windowed_mean = window_mean(RMS_windowed,mean_window);
RMS_robust_mean = window_mean(RMS_robust,mean_window);

figure
subplot(2,1,1)
plot(n(1:N*P/2),abs(U(1:N*P/2)),'.','Markersize',15,'Linewidth',2)
xlabel("DFT bins (k)")
title("|U| (linear)")
plot_options(gca)
subplot(2,1,2)
plot(n(1:N*P/2),db(Y(1:N*P/2)),'.','Markersize',15,'Linewidth',2)
xlabel("DFT bins (k)")
title("|Y|_{dB}")
plot_options(gca)


figure
hold on
plot(fexc/Fs,db(G_actual_exc/10),'k-','Linewidth',2,'Displayname',"|G/10|")
% plot(fexc/Fs,db(RMS_normal),'r.--','Markersize',15,'Displayname',"RMS[|G_{est}-G|] (No preprocessing)")
% plot(fexc/Fs,db(RMS_robust),'g.:','Markersize',15,'Displayname',"RMS[|G_{est}-G|] (Robust LPM)")
% plot(fexc/Fs,db(RMS_windowed),'mo:','Markersize',8,'Displayname',"RMS[|G_{est}-G|] (Hann window)")
plot(fexc_mean/Fs,db(RMS_normal_mean),'rx--','Markersize',8,'Displayname',"RMS(|G_{est}-G|) (No preprocessing)")
plot(fexc_mean/Fs,db(RMS_robust_mean),'b.:','Markersize',15,'Displayname',"RMS(|G_{est}-G|) (Robust LPM)")
plot(fexc_mean/Fs,db(RMS_windowed_mean),'mo:','Markersize',6,'Displayname',"RMS(|G_{est}-G|) (Hann window)")
xlabel('Normalized frequency')
plot_options(gca)
title("Magnitude (dB)")
ylim([-20,5])
legend('Location',[0.4925    0.7123    0.4786    0.1905])
if NR_trials == 100
    print(gcf,'figures/LPM_example','-depsc')
end



figure
% plot(fexc/Fs,diffRMS,'b:.','Markersize',15)
plot(fexc_mean/Fs,diffRMS_mean,'b:.','Markersize',15)
xlabel('Normalized frequency')
plot_options(gca)
ylabel("Magnitude (dB)")
title("RMS Robust LPM (dB) - RMS Hann window (dB)")


figure
hold on
plot(fall/Fs,db(Gall),'k-','Linewidth',2,'Displayname',"|G|")
plot(fexc/Fs,db(G_est_robust(:,1)),'b.-','Markersize',15,'Displayname',"|G_{est}| (Robust LPM)")
plot(fexc/Fs,db(varG(:,1)),'rx--','Markersize',15,'Displayname',"|Var(G_{est})| (Robust LPM)")
xlabel('Normalized frequency')
plot_options(gca)
title("Magnitude (dB)")
legend('Location','Best')

function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end

function x_mean = window_mean(x,windowsize)
    N = numel(x);
    windows = floor(N/windowsize);
    for i = windows:-1:1
        x_mean(i) = mean(x((1:windowsize)+(i-1)*windowsize));
    end
end