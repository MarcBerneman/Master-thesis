classdef TD_controller
    %TD_CONTROLLER Find optimal controller using TD methods
    
    properties
        np % scalar
        M % TF
        F % TF or a scalar
        beta % TF
        Phir_exact % TF or scalar
    end
    
    methods
        function obj = TD_controller(np,M,F,beta,Phir_exact)
            %TD_CONTROLLER Construct an instance of this class
            obj.np = np;
            obj.M = M;
            obj.F = F;
            obj.beta = beta ;
            obj.Phir_exact = Phir_exact;
        end
        
        function [rho,Delta] = optimize(obj,u,y,l1,deltaN)
            % r: N x 1 vector
            % y: N x 1 vector
            % l1: scalar
            % delta_N: scalar smaller than 1 or [] for unconstrained optimization
            NP = numel(u);
            N = NP/obj.np;
            if even(N)
                k_max = N/2-1;
            else
                k_max = (N-1)/2;
            end

            W = obj.F*(1-obj.M)/obj.Phir_exact;
            uW = my_dlsim(W,u);

            % epsilon calculation
            n_rho = size(obj.beta,1)-1;
            rho = sdpvar(n_rho+1,1,'full');
            epsilon = my_dlsim(obj.M,u) - my_dlsim((1-obj.M)*obj.beta,y)*rho;

            % cost function
            J = cost(obj,uW,epsilon,l1);
            Rr = calc_correlation_TD(u,u,0:N-1);
            Rrepsilon = calc_correlation_TD(u,epsilon,0:N-1);
            Delta = fft(Rrepsilon)./fft(Rr);
            if exist('deltaN','var')
                constraints = abs(Delta(2:k_max+1)) <= deltaN; %No DC
            else
                constraints = [];
            end

            %optimize
            evalc("optimize(constraints,J);");
            rho = value(rho);
            Delta = value(Delta);
        end
        
        function J = cost(obj,uW,epsilon,l1)
            corr_func = calc_correlation_TD(uW,epsilon,-l1:l1);
            J = corr_func'*corr_func;
        end
        
        function [rho,Delta] = fast_optimize(obj,u,y,l1) 
            W = obj.F*(1-obj.M)/obj.Phir_exact;
            uW = my_dlsim(W,u);

            % epsilon calculation
            uM = my_dlsim(obj.M,u);
            y2 = my_dlsim((1-obj.M)*obj.beta,y);
            lags = -l1:l1;
            a = calc_correlation_TD(uW,uM,lags);
            for i = size(y2,2):-1:1
                D(:,i) = calc_correlation_TD(uW,y2(:,i),lags);
            end
            rho = D\a;
   
            Delta = 0;
        end
    end
end

