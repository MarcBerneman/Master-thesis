function [G_est,varG] = nonparametric_estimate(r0,u,y,P,Ts,ExcitedHarm,order,dof,estimate_transient)
    %NONPARAMETRIC_ESTIMATE Estimate the TF nonparametrically by using
    % the Robust Local Polynomial Method. 
    % y is the measured output of the plant: NP x 1 vector
    % u is the measured input of the plant: NP x 1 vector
    % r0 is the reference signal: NP x 1 vector
    % Ts is the sampling period: scalar
    % P is the number of periods in u and y: scalar
    % ExcitedHarm are the DFT bins that are excited in the input: K x 1 vector
    % order is the order used for the Robust Local Polynomial Method: scalar
    % dof is the degrees of freedom used for the Robust Local Polynomial Method: scalar
    % estimate_transient: boolean
    
    N = numel(u)/P;
    data.r(1,1,1,:) = r0;
    data.u(1,1,1,:) = u;
    data.y(1,1,1,:) = y;
    data.Ts = Ts;
    data.N = N;
    data.ExcitedHarm = ExcitedHarm;
    
    if estimate_transient
        method.transient = 1;
    else
        method.transient = 0;
    end
    method.dof = dof;
    method.order = order;
    [~, Z, ~, G_est, CvecG, ~, ~] = RobustLocalPolyAnal(data, method);
    Z = squeeze(Z);
    G_est = squeeze(G_est);
    varG = squeeze(CvecG.n);
end

