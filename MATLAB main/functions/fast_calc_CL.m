function CLk = fast_calc_CL(Gk,betak,rho)
    for i = size(rho,2):-1:1
        Kk = betak.' * rho(:,i);
        CLk(:,i) = Kk.*Gk./(1+Kk.*Gk);
    end
end

