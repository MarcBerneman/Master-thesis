clear
close all

Ts = 0.05;
Fs = 1/Ts;
num = [0,0,0,0.7893];
den = [1,-1.418,1.59,-1.316,0.886];
G = tf(num,den,Ts,'Variable','q^-1');
qinv = tf([0 1],1,Ts,'Variable','q^-1');
alpha = 0.606;
M = qinv^3*(1-alpha)^2/(1-alpha*qinv)^2;

fprintf("DC gain G = %2.2f\n",dcgain(G))
fprintf("DC gain M = %2.2f\n",dcgain(M))

N = 1000;
if even(N)
    k_max = N/2-1;
else
    k_max = (N-1)/2;
end
f = (0:k_max)*Fs/N;
Gbode = squeeze(freqresp(G,f,'Hz'));
figure
plot(f,db(Gbode))
legend("G")

Mbode = squeeze(freqresp(M,f,'Hz'));
figure
plot(f,db(Mbode))
legend("M")

save('System_data','M','G','Ts','qinv')