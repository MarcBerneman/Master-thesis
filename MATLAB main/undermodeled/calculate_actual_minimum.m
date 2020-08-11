clear
close all

N = 255; %samples per period
F_all = 1;
NR_harmonics = 127;
ExcitedHarm = 1:NR_harmonics;

load('System_data')
Fs = 1/Ts;

% Controller
n_rho = 4;
qinv = tf([0,1],1,Ts,"variable","q^-1");
for i = n_rho:-1:0
    beta_prime(i+1,1) = qinv^i; %better for some computations
    beta(i+1,1) = qinv^i / (1-qinv);
end

k_all = (0:N-1).';
f_all = k_all*Fs/N;
% frequence response for perfect filtering
M_all = squeeze(freqresp(M,f_all,'Hz'));
beta_all = squeeze(freqresp(beta,f_all,'Hz'));
if n_rho > 0
    beta_all = beta_all.';
end

FD_method = FD_controller(M_all,F_all,beta_all);
fk = f_all(ExcitedHarm + 1);
Gk = squeeze(freqresp(G,fk,'Hz'));
[rho_actual,~] = FD_method.fast_optimize_no_l1(Gk,ExcitedHarm);

save('rho_actual','rho_actual')
