function [Jmr,J] = fast_calc_cost(Mk,Gk,betak,rho)
   for i = size(rho,2):-1:1
        Kk = betak.' * rho(:,i);
        CLk = Kk.*Gk./(1+Kk.*Gk);

        bode_mr = Mk-CLk;
        bode = (1-Mk).*(Mk-Kk.*(1-Mk).*Gk);

        Jmr(i) = bode_mr'*bode_mr;
        J(i) = bode'*bode;
    end
    Jmr = Jmr/numel(bode_mr);
    J = J/numel(bode);
end

