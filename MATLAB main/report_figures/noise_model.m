clear
close all

%% load data

load("../Noise_filters")

N = 255;
if even(N)
    k_max = N/2-1;
else
    k_max = (N-1)/2;
end
ExcitedHarm = 1:k_max;

Ts = 1;
Fs = 1/Ts;
f = (0:k_max)*Fs/N;
f_exc = f(ExcitedHarm+1);
Sy = squeeze(freqresp(Sy,f_exc,'Hz'));
Sy_flat = ones(size(Sy));
%% figures
figure
hold on
plot(f_exc/Fs,db(Sy),'b-','Linewidth',4,'Markersize',15,'Displayname','Coloured noise')
plot(f_exc/Fs,db(Sy_flat),'r:','Linewidth',4,'Markersize',15,'Displayname','White noise')
legend("Location","Best")
grid on
xlabel("Normalized frequency")
title("Magnitude (dB)")
plot_options(gca)
print(gcf,"figures/noise_models","-depsc")

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end