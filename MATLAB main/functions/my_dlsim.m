function y = my_dlsim(G,u)
    % Simulate linear discrete time system with zero initial conditions. 
    % Input u can be an sdvpar.
    % G: SIMO (multiple outputs possible)
    % u: N x 1
    N = numel(u);
    Ts = G.Ts;
    g = impulse(G,N*Ts);
    y = Ts*conv2(u,g);
    y = y(1:N,:);
end

