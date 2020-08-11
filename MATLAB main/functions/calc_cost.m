function [Jmr,J] = calc_cost(M,G,beta,rho,N,ExcitedHarm)
    if even(N)
        k_max = N/2-1;
    else
        k_max = (N-1)/2;
    end
    Fs = 1/M.Ts;
    f = (0:k_max)*Fs/N;
    f_exc = f(ExcitedHarm+1);
    
    for i = size(rho,2):-1:1
        K = beta.' * rho(:,i);
        CL = feedback(K*G,1);

        bode_mr = squeeze(freqresp(M-CL,f_exc,'Hz'));
        bode = squeeze(freqresp((1-M)*(M-K*(1-M)*G),f_exc,'Hz'));

        Jmr(i) = bode_mr'*bode_mr;
        J(i) = bode'*bode;
    end
end

