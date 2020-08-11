function corr_func = calc_correlation_TD(a,b,lags)
    row = numel(lags); %start at end to create matrix with correct size
    N = numel(a);
    for tau = lags
        lower_t = max([1,tau+1]);
        upper_t = min([N,N+tau]);
        t = lower_t:upper_t;
        corr_func(row,1) = 1/N*sum(a(t-tau).*b(t));
        row = row - 1;
    end
    corr_func = flipud(corr_func);
end

