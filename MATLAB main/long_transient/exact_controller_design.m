clear
close all

Ts = 2;
Fs = 1/Ts;

w0 = 2*pi*0.3;
zeta = 0.01;

s = tf('s');
Gs = w0^2/(s^2 + 2*zeta*w0*s + w0^2);

Gz = c2d(Gs,Ts,'zoh');

figure
bodemag(Gz,{2*pi*Ts/255,pi*Ts})

figure
impulse(Gz)

z = tf('z',Ts);
Kp = 0.1;
Ki = 0.1;
K_exact = Kp+Ki*Ts/(z-1);

rho_actual = tfdata(K_exact);
rho_actual = rho_actual{1}.';

M_exact = feedback(K_exact*Gz,1);

fprintf("DC gain G = %2.2f\n",dcgain(Gs))
fprintf("DC gain M = %2.2f\n",dcgain(M_exact))

figure
bodemag(M_exact,{2*pi*Ts/255,pi*Ts})
legend("M_{exact}")

figure
hold on
step(Gz)
step(M_exact)
legend("G","M_{exact}")

G = Gz;
save('System_data','M_exact','K_exact','G','Ts','z')
save('rho_actual','rho_actual')