clear
close all

Ts = 1;
Fs = 1/Ts;

qinv = tf([0 1],1,Ts,'Variable','q^-1');
a = -0.75;
G = tf([0,1+a],[1,a],Ts,'Variable','q^-1');

figure
bodemag(G,{2*pi*Ts/255,pi*Ts})

Kp = 0.6;
Ki = 0.3;
K_exact = Kp+Ki*qinv/(1-qinv);

rho_actual = tfdata(K_exact);
rho_actual = rho_actual{1}.';

M_exact = feedback(K_exact*G,1);

fprintf("DC gain G = %2.2f\n",dcgain(G))
fprintf("DC gain M = %2.2f\n",dcgain(M_exact))

figure
bodemag(M_exact,{2*pi*Ts/255,pi*Ts})
legend("M_{exact}")

figure
step(M_exact)
legend("M_{exact}")

save('System_data','M_exact','K_exact','G','Ts','qinv')
save('rho_actual','rho_actual')