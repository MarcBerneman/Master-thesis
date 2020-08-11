clear
close all

load("TF2_6")
load("K")
load("M")

CL = feedback(TF*K,1);

Fs = 78.125e3;
N = 2048;
f = (0:N-1)*Fs/N;
fmax = Fs/2;
kmax = floor(fmax/Fs*N);
ExcitedHarm = 1:kmax;

f = f(ExcitedHarm+1);
CL_bode = squeeze(freqresp(CL,f,'Hz'));
TF_bode = squeeze(freqresp(TF,f,'Hz'));
%% Figures
figure
hold on
plot(f/1000,db(TF_bode),'b-','linewidth',2,'Displayname',"G_{BLA}")
plot(f/1000,db(CL_bode),'g-.','linewidth',2,'Displayname',"CL")
xlabel('f (kHz)')
ylabel('Magnitude (dB)')
set(gca,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',12)
legend('Location','Best')
grid on

print("real_system_parametric_closed_loop",'-djpeg')

figure
pzmap(CL)
print("real_system_parametric_poles_zeros",'-djpeg')

figure
hold on
step(CL,(0:200)*0.00001)
step(M,(0:200)*0.00001)