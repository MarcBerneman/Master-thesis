clear
close all

%% load data
for study_case = ["simple","long_transient","undermodeled"]
    clearvars -except study_case
    load("../" + study_case + "/System_data")
    if ~exist("M","Var")
        M = M_exact;
    end

    N = 255;
    if even(N)
        k_max = N/2-1;
    else
        k_max = (N-1)/2;
    end
    ExcitedHarm = 1:k_max;
    Fs = 1/M.Ts;
    f = (0:k_max)*Fs/N;
    f_exc = f(ExcitedHarm+1);
    Gk = squeeze(freqresp(G,f_exc,'Hz'));
    Mk = squeeze(freqresp(M,f_exc,'Hz'));

    %% figures
    figure
    hold on
    plot(f_exc/Fs,db(Gk),'b-','Linewidth',4,'Markersize',15,'Displayname','G')
    plot(f_exc/Fs,db(Mk),'r:','Linewidth',4,'Markersize',15,'Displayname','M')
    legend("Location","Best")
    grid on
    xlabel("Normalized frequency")
    title("Magnitude (dB)")
    plot_options(gca)
    print(gcf,"figures/G_and_M_"+study_case,"-depsc")
end

%%
load("../simple/System_data")
G_simple = G;
Ts_simple = G_simple.Ts;
load("../long_transient/System_data")
G_long = G;
Ts_long = G_long.Ts;

[g_long,t] = impulse(G_long);
n = t/G_long.Ts;
g_simple = impulse(G_simple,n*Ts_simple);

figure
hold on
s2 = stairs(n,g_long,'r-','Linewidth',2,'Markersize',15,'Displayname','Long transient system');
s1 = stairs(n,g_simple,'b:','Linewidth',4,'Markersize',15,'Displayname','Simple system');
legend([s1,s2],"Location","Best")
grid on
xlabel("n (samples)")
title("Impulse response")
xlim([-Inf,60])
plot_options(gca)
print(gcf,"figures/impulse_response_simple_long","-depsc")

function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',18)
    grid on
end