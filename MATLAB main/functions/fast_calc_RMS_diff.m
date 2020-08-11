function RMS = fast_calc_RMS_diff(Mk,Gk,betak,rho)
    RMS = zeros(numel(Gk),1);
    for i = size(rho,2):-1:1
        Kk = betak.' * rho(:,i);
        CLk = Kk.*Gk./(1+Kk.*Gk);

        bode_mr = Mk-CLk;
        RMS = RMS + abs(bode_mr).^2;
    end
    RMS = sqrt(RMS/size(rho,2));
end

