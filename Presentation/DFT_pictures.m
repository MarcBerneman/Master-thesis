clear
close all

fs = 1;
ts = 1/fs;
N = 100;
t_sampled = (0:N-1)*ts;
OS = 10;
t_cont = (0:N*OS-1)*(ts/OS);

for f1 = [0.1,0.4,0.6]
    x = @(t) sin(2*pi*f1*t);
    x_cont = x(t_cont);
    x_sampled = x(t_sampled); 

    X_cont = 1/(N*OS)*fft(x_cont);
    X_sampled = 1/N*fft(x_sampled);
    f_cont = (0:N*OS-1)*(fs*OS)/(N*OS);
    f_sampled = (0:N-1)*fs/N;

    X_cont = X_cont(1:end/2);
    X_sampled = X_sampled(1:end/2);
    f_cont = f_cont(1:end/2);
    f_sampled = f_sampled(1:end/2);

    figure
    hold on
    plot(t_cont,x_cont,'b-','Linewidth',2)
    plot(t_sampled,x_sampled,'r.','Markersize',15)
    grid on

    figure('Position',[258,246,705,420])
    hold on
    stem(f_cont(abs(X_cont) > 1e-8),abs(X_cont(abs(X_cont) > 1e-8)),'bo','Linewidth',1.5,'Markersize',7,'Displayname','Actual')
    stem(f_sampled(abs(X_sampled) > 1e-8),abs(X_sampled(abs(X_sampled) > 1e-8)),'r.','Linewidth',1,'Markersize',13,'Displayname','Sampled')
    xline(fs/2,'Linewidth',2,'Displayname','Nyquist frequency');
    grid on
    xlim([0,fs])
    legend('Location','Northeast','Fontsize',12)
    ylim([0,0.6])
    % xlabel("")
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlabel('Frequency')
    ylabel('Magnitude')
    set(gca,'Fontsize',12)
    print(gcf,"pictures/nyquist_shannon_" + f1 + ".jpg","-djpeg")
end

%%
f1 = 0.05;
x = @(t) sin(2*pi*f1*t);

for N1 = [20,25]
t1 = (0:N1-1)*ts;
f1 = (0:N1-1)*fs/N1;
x1 = x(t1);
X1 = 1/N1*fft(x1);

figure
subplot(2,1,1)
plot(t1,x1,'r.--','Markersize',15)
xlabel('Time')
ylabel('Signal')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'Fontsize',12)
subplot(2,1,2)
stem(f1,abs(X1),'r.','Markersize',13)
xlim([0,0.5])
xlabel('Frequency')
ylabel('Magnitude')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'Fontsize',12)
print(gcf,"pictures/leakage_" + N1 + ".jpg","-djpeg")
end