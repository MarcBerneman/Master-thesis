classdef FD_WNLS_controller
    %FD_WNLS_CONTROLLER Find optimal controller using FD methods and a
    % weighted nonlinear least squares cost function. 
    % Non-parametric transfer function determined using the Robust Local 
    % Polynomial method (order = 2)
    
    properties
        M_all % T x 1 vector
        F_all % T x 1 vector or a scalar
        beta_all % T x (nrho+1) matrix
    end
    
    methods
        function obj = FD_WNLS_controller(M_all,F_all,beta_all)
            %FD_CONTROLLER Construct an instance of this class
            obj.M_all = M_all;
            obj.F_all = F_all;
            obj.beta_all = beta_all;
        end
        
        function [rho,Delta] = optimize(obj,G_est,varG,ExcitedHarm,initEst,maxiter,precision)
            % G_est: K x 1 vector
            % varG: K x 1 vector
            % ExcitedHarm: K x 1 vector
            % initEst: nrho x 1 vector or "LS"
            % maxiter: scalar
            % precision: scalar
            a = obj.F_all.*(1-obj.M_all).*obj.M_all;
            B = obj.F_all.*(1-obj.M_all).^2.*obj.beta_all;
            a = a(ExcitedHarm+1);
            B = B(ExcitedHarm+1,:);
            
            if strcmp(initEst,"LS")
                rho = LS_solution(obj,G_est,ExcitedHarm);
            else
                rho = initEst;
            end
            for i = 1:maxiter-1
                dJ = WNLS_cost_gradient(obj,a,B,G_est,varG,rho);
                dJ2 = calc_WNLS_cost_hessian(obj,a,B,G_est,varG,rho);
                rho_next = rho - dJ2\dJ;
                diff = abs((rho_next-rho)./rho);
                rho = rho_next;
                if all(diff < precision)
                    break
                end
            end
            
            Delta = obj.M_all(ExcitedHarm+1) - ...
                (1-obj.M_all(ExcitedHarm+1)).*G_est.*obj.beta_all(ExcitedHarm+1,:)*rho;
        end
        
        function rho = LS_solution(obj,G_est,ExcitedHarm)
            a = obj.F_all.*(1-obj.M_all).*obj.M_all;
            B = obj.F_all.*(1-obj.M_all).^2.*obj.beta_all;
            a = a(ExcitedHarm+1);
            B = B(ExcitedHarm+1,:);
            right = zeros(size(B,2),1);
            left = zeros(size(B,2));
            for k = 1:numel(G_est)
                right = right + B(k,:)'*conj(G_est(k))*a(k);
                left = left + B(k,:)'*B(k,:)*abs(G_est(k))^2;
            end
            rho = real(left)\real(right);
        end
              
        function J = WNLS_cost(obj,a,B,G_est,varG,rho)
            b = B*rho;
            b2 = abs(b).^2;
            H2 = abs(a-b.*G_est).^2;
            varH = varG.*b2;
            J = sum(H2./varH);
        end

        function dJ = WNLS_cost_gradient(obj,a,B,G_est,varG,rho)
            b = B*rho;
            dJ = zeros(size(rho));
            for k = 1:numel(G_est)
                dJ = dJ + 2*real(B(k,:)'*b(k)*conj(a(k))*(b(k)*G_est(k)-a(k)))/(varG(k)*abs(b(k))^4);
            end
        end
        
        function dJ2 = calc_WNLS_cost_hessian(obj,a,B,G_est,varG,rho)
            b = B*rho;
            dJ2 = zeros(numel(rho));
            denom = abs(b).^4;
            for k = 1:numel(G_est)
                d_denom = 4*abs(b(k))^2*real(B(k,:)'*b(k));
                for i = 1:numel(rho)
                    Bi = B(k,i);
                    num = 2*real(conj(Bi)*b(k)*conj(a(k))*(b(k)*G_est(k)-a(k)));
                    d_num = 4*real(conj(Bi*a(k))*G_est(k)*B(k,:).'*b(k)) - ...
                            2*real(conj(Bi)*abs(a(k))^2*B(k,:).');
                    dJ2(:,i) = dJ2(:,i) + (d_num.*denom(k)-num.*d_denom)/(varG(k)*denom(k)^2);
                end
            end
        end        
        

    end
end

