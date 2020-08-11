clear
close all
rng(1)

N = 40;
n = 0:N-1;
Ts = 0.5; Fs = 1/Ts;
fmin = Fs/N;
fmax = 3*Fs/N;
df = Fs/N;

RMS = 1/8.165;
[x,k_exc] = generate_MS(N,Fs,RMS,fmin,fmax,df);
X = 1/N*fft(x);

P = 3;
n_rep = 0:N*P-1;
x_rep = repmat(x,P,1);
X_rep = 1/(N*P)*fft(x_rep);

figure('Units','Normalized','Position',[0.1,0.1,0.6,0.8])
subplot(2,3,1)
plot(n,x,'-','Linewidth',2)
xlabel('Samples (n)')
title("x(n)")
plot_options(gca)
subplot(2,3,4)
stem(n(1:N/2),abs(X(1:N/2)),'Linewidth',2)
xlabel('DFT bins (k)')
title("|X|")
plot_options(gca)
subplot(2,3,[2,3])
plot(n_rep,x_rep,'-','Linewidth',2)
xlabel('Samples (n)')
title("x(n) (" + P + " periods)")
plot_options(gca)
subplot(2,3,[5,6])
stem(n_rep(1:N*P/2),abs(X_rep(1:N*P/2)),'Linewidth',2)
xlabel('DFT bins (k)')
title("|X| (" + P + " periods)")
plot_options(gca)
print(gcf,'figures/periodic','-depsc')

%% output
sigma_y = 0.005;

w0 = 2*pi*0.3;
zeta = 0.01;
s = tf('s');
G = w0^2/(s^2 + 2*zeta*w0*s + w0^2);
z = tf('z',Ts);
Gz = c2d(G,Ts,'zoh');


ny = sigma_y*randn(size(n_rep));
y_rep0 = lsim(Gz,x_rep,n_rep*Ts).';
y_rep = y_rep0 + ny;
Y_rep0 = 1/(N*P)*fft(y_rep0);
Y_rep = 1/(N*P)*fft(y_rep);

figure('Units','Normalized','Position',[0.1,0.1,0.6,0.8])
subplot(2,1,1)
stem(n_rep(2:N*P/2),abs(X_rep(2:N*P/2)),'Linewidth',2)
xlabel('DFT bins (k)')
title("|U|")
plot_options(gca)
subplot(2,1,2)
hold on
plot(n_rep(2:N*P/2),db(Y_rep0(2:N*P/2)),'.','Markersize',15,'Linewidth',2,'Displayname','Without noise')
plot(n_rep(2:N*P/2),db(Y_rep(2:N*P/2)),'o','Markersize',8,'Linewidth',2,'Displayname',"With noise (\sigma_y = " + sigma_y + ")")
xlabel('DFT bins (k)')
title("|Y|_{dB}")
legend
annotation('textarrow',[0.15284046358717 0.20284046358717],...
    [0.497919550760044 0.437919550760044],'TextEdgeColor',[0 0 0],...
    'String','G(kP) U(kP)',...
    'LineWidth',1,...
    'FontSize',12);
annotation('textarrow',[0.444444444444444 0.366300366300366],...
    [0.403908794788274 0.379478827361564],'TextEdgeColor',[0 0 0],...
    'String','T(kP+r)',...
    'LineWidth',1,...
    'FontSize',12);
annotation('textarrow',[0.577533577533577 0.563255192830595],...
    [0.379478827361564 0.28679445745114],'TextEdgeColor',[0 0 0],...
    'String','N_Y(kP+r)',...
    'LineWidth',1,...
    'FontSize',12);
plot_options(gca)
print(gcf,'figures/periodic_output','-depsc')
%% functions
function plot_options(gca)
    set(gca,'Linewidth',2)
    set(gca,'Fontsize',12)
    grid on
end
