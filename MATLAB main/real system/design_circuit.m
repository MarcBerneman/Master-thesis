clear
close all

load("rho")
P.gain = rho(1);
I.gain = rho(2);
D1.gain = rho(3);
D2.gain = rho(4)/rho(3);

%% P
P.R1 = 1e3;
P.R2 = P.R1 * P.gain;

%% I
I.f0 = 10;
I.R1 = 220;
I.C = 1/(I.R1*I.gain);
I.R2 = 1/(2*pi*I.C*I.f0);
I.DCGAIN_dB = db(I.R2/I.R1);
%% D1
% D1.R = 1e3;
% D1.C = D1.gain/D1.R;
D1.C = 150e-9;
D1.R = D1.gain/D1.C;

%% D2
% D2.R = 1e3;
% D2.C = D2.gain/D2.R;
% D2.C = 120e-9;
% D2.R = D2.gain/D2.C;
D2.C = 120e-9;
D2.Rf = D2.gain/D2.C;
D2.maxgain = 10;
D2.Rin = D2.Rf/D2.maxgain;
D2.cutoff_point = 100e3;
D2.Cf = D2.cutoff_point*2*pi/D2.Rin;

load("K")

f = logspace(log10(1),log10(1e6));
Kk = squeeze(freqresp(K,f,'Hz'));
Pk = ones(size(Kk))*P.gain;
D1k = squeeze(freqresp(D1.gain*tf('s'),f,'Hz'));
D2k = squeeze(freqresp(D2.gain*tf('s'),f,'Hz'));
Ik = squeeze(freqresp(I.gain/tf('s'),f,'Hz'));
DDk = D1k.*D2k;

file = "Circuit files/Vout_with_subtractor_real_values_more_freqs.txt";
[Ksim,fsim] = load_sim_data(file);
file = "Circuit files/Vout_with_subtractor_real_values_more_freqs_no_finite_D2.txt";
[Ksim_no_finite_D2,~] = load_sim_data(file);
% [D1sim,DDsim,Isim,Psim,fsim] = load_sim_data2();

figure
hold on
plot(f,db(Kk),'k-','linewidth',4,'Displayname',"K")
plot(fsim,db(Ksim),'b-.','linewidth',4,'Displayname',"K_{sim}")
plot(fsim,db(Ksim_no_finite_D2),'r:','linewidth',3,'Displayname',"K_{sim} (no added components)")

% plot(f,db(Pk+D1k+DDk+Ik),'linewidth',2,'Displayname',"Sum")
% plot(fsim,db(Psim+D1sim+DDsim+Isim),'linewidth',2,'Displayname',"Sum_{sim}")

% plot(f,db(D1k),'linewidth',2,'Displayname',"D1")
% plot(fsim,db(D1sim),'--','linewidth',2,'Displayname',"D1_{sim}")

% plot(f,db(D2k),'linewidth',2,'Displayname',"D2")

% plot(f,db(Pk),'linewidth',2,'Displayname',"P")
% plot(fsim,db(Psim),'--','linewidth',2,'Displayname',"P_{sim}")

% plot(f,db(Ik),'linewidth',2,'Displayname',"I")
% plot(fsim,db(Isim),'--','linewidth',2,'Displayname',"I_{sim}")

% plot(f,db(DDk),'g:','linewidth',2,'Displayname',"D1.D2")
% plot(fsim,db(DDsim),'--','linewidth',2,'Displayname',"D1.D2_{sim}")
xlabel('f (Hz)')
ylabel('Magnitude (dB)')
set(gca,'linewidth',2)
set(gca,'xscale','log')
xlim([min(f),max(f)])
xticks(logspace(0,6,7))
plot_options(gca)
grid off
legend('Location','Best','fontsize',14)
title("Optimal controller vs. simulated controller")
print(gcf,"controller_simulation",'-depsc')

% figure
% hold on
% plot(f,angle(Kk)*180/pi,'k-','linewidth',2,'Displayname',"K")
% plot(fsim,angle(Ksim)*180/pi,'b-.','linewidth',2,'Displayname',"K_{sim}")

% plot(f,angle(Pk+D1k+DDk+Ik)*180/pi,'linewidth',2,'Displayname',"Sum")
% plot(fsim,angle(Psim+D1sim+DDsim+Isim)*180/pi,'linewidth',2,'Displayname',"Sum_{sim}")

% plot(f,angle(D1k)*180/pi,'linewidth',2,'Displayname',"D1")
% plot(fsim,angle(D1sim)*180/pi,'--','linewidth',2,'Displayname',"D1_{sim}")

% plot(f,angle(D2k)*180/pi,'linewidth',2,'Displayname',"D2")

% plot(f,angle(Pk)*180/pi,'linewidth',2,'Displayname',"P")
% plot(fsim,angle(Psim)*180/pi,'--','linewidth',2,'Displayname',"P_{sim}")

% plot(f,angle(Ik)*180/pi,'linewidth',2,'Displayname',"I")
% plot(fsim,angle(Isim)*180/pi,'--','linewidth',2,'Displayname',"I_{sim}")

% plot(f,angle(DDk)*180/pi,'g:','linewidth',2,'Displayname',"D1.D2")
% plot(fsim,angle(DDsim)*180/pi,'--','linewidth',2,'Displayname',"D1.D2_{sim}")
% xlabel('f (Hz)')
% ylabel('Angle (°)')
% set(gca,'linewidth',2)
% set(gca,'xscale','log')
% xlim([min(f),max(f)])
% set(gca,'linewidth',2)
% set(gca,'fontsize',12)
% legend('Location','Best')
% grid on

function [Ksim,fsim] = load_sim_data(file)
    data = readmatrix(file,'NumHeaderLines',1,'DecimalSeparator','.','Delimiter',{'\t',','});
    fsim = data(:,1);
    Ksim = data(:,2) + 1j*data(:,3);
end

function [D1sim,DDsim,Isim,Psim,fsim] = load_sim_data2()
    file = "Circuit files/PIDD.txt";
    data = readmatrix(file,'NumHeaderLines',1,'DecimalSeparator','.',...
        'Delimiter',{'\t',',','\t',',','\t',',','\t',','});
    fsim = data(:,1);
    D1sim = data(:,2) + 1j*data(:,3);
    DDsim = data(:,4) + 1j*data(:,5);
    Isim = data(:,6) + 1j*data(:,7);
    Psim = data(:,8) + 1j*data(:,9);
end
    
function plot_options(gca)
    set(gca,'Linewidth',4)
    set(gca,'Fontsize',14)
    grid on
end