classdef FD_controller
    %FD_CONTROLLER Find optimal controller using FD methods
    % Non-parametric transfer function determined using the Robust Local 
    % Polynomial method (order = 2)
    
    properties
        M_all % T x 1 vector
        F_all % T x 1 vector or a scalar
        beta_all % T x (nrho+1) matrix
    end
    
    methods
        function obj = FD_controller(M_all,F_all,beta_all)
            %FD_CONTROLLER Construct an instance of this class
            obj.M_all = M_all;
            obj.F_all = F_all;
            obj.beta_all = beta_all;
        end
        
        function [rho,Delta] = optimize(obj,G_est,ExcitedHarm,l1,deltaN)
            % G_est: K x 1 vector
            % ExcitedHarm: K x 1 vector
            % l1: scalar
            % delta_N: scalar smaller than 1 or [] for unconstrained optimization
            N = numel(obj.M_all);
            if even(N)
                k_max = N/2-1;
            else
                k_max = (N-1)/2;
            end
            G_est_ext = zeros(k_max,1);
            G_est_ext(ExcitedHarm) = G_est;
            if even(N)
                G_est_ext = [eps;G_est_ext;eps;flipud(conj(G_est_ext))];
            else
                G_est_ext = [eps;G_est_ext;flipud(conj(G_est_ext))];
            end
            
            % optimization variables
            n_rho = size(obj.beta_all,2)-1;
            rho = sdpvar(n_rho+1,1,'full');
            Delta = obj.M_all - (1-obj.M_all).*G_est_ext.*obj.beta_all*rho;
            
            % cost function
            J = cost(obj,G_est_ext,obj.M_all,Delta,obj.F_all,l1);

            if exist('deltaN','var')
                constraints = abs(Delta(ExcitedHarm+1)) <= deltaN; %No DC
            else
                constraints = [];
            end

            %optimize
            evalc("optimize(constraints,J);");
            rho = value(rho);
            Delta = value(Delta);
        end
                
        function J = cost(obj,G_est_ext,M_all,Delta,F_all,l1)
            N = numel(G_est_ext);
            if even(N)
                k_max = N/2-1;
            else
                k_max = (N-1)/2;
            end
            H = F_all.*(1-M_all).*Delta;
            H(1) = 0; %make zero mean
            if even(N)
                H(k_max+1) = 0;
            end
%             h = fftshift(ifft(H.*conj(H)));
%             corr_func = h((l1:-1:-l1)+k_max+1);
%             J = abs(corr_func'*corr_func);
            h = ifft(H);
            J = h(1:l1+1)'*h(1:l1+1);
%             J = h(1:2*l1+1)'*h(1:2*l1+1);
        end
        
        function [rho,Delta] = fast_optimize(obj,G_est,ExcitedHarm,l1)
            % G_est: K x 1 vector
            % ExcitedHarm: K x 1 vector
            % l1: scalar
            N = numel(obj.M_all);
            if even(N)
                k_max = N/2-1;
            else
                k_max = (N-1)/2;
            end
            
            a = obj.F_all.*(1-obj.M_all).*obj.M_all;
            B = obj.F_all.*(1-obj.M_all).^2.*obj.beta_all;
            a = a(ExcitedHarm+1);
            B = B(ExcitedHarm+1,:);
            D = B.*G_est;
            
            a_ext = zeros(k_max,1);
            a_ext(ExcitedHarm) = a;
            if even(N)
                a_ext = [0;a_ext;0;flipud(conj(a_ext))];
            else
                a_ext = [0;a_ext;flipud(conj(a_ext))];
            end
            
            D_ext = zeros(k_max,size(D,2));
            D_ext(ExcitedHarm,:) = D;
            if even(N)
                D_ext = [zeros(1,size(D,2));D_ext;zeros(1,size(D,2));flipud(conj(D_ext))];
            else
                D_ext = [zeros(1,size(D,2));D_ext;flipud(conj(D_ext))];
            end
            
            aTD = ifft(a_ext);
            DTD = ifft(D_ext);
            rho = DTD(1:l1+1,:)\aTD(1:l1+1);
            Delta = 0;
        end
        
        function [rho,Delta] = fast_optimize_no_l1(obj,G_est,ExcitedHarm)
            a = obj.F_all.*(1-obj.M_all).*obj.M_all;
            B = obj.F_all.*(1-obj.M_all).^2.*obj.beta_all;
            a = a(ExcitedHarm+1);
            B = B(ExcitedHarm+1,:);
            D = B.*G_est;
            aRI = [real(a);imag(a)];
            DRI = [real(D);imag(D)];
            rho = DRI\aRI;
            
            Delta = obj.M_all(ExcitedHarm+1) - (1-obj.M_all(ExcitedHarm+1)).*G_est.*obj.beta_all(ExcitedHarm+1,:)*rho;
        end
        
        function [rho,Delta] = optimize_no_l1(obj,G_est,ExcitedHarm,deltaN)
            % optimization variables
            n_rho = size(obj.beta_all,2)-1;
            rho = sdpvar(n_rho+1,1,'full');
            Mk = obj.M_all(ExcitedHarm+1);
            betak = obj.beta_all(ExcitedHarm+1,:);
            
            Delta = Mk - (1-Mk).*G_est.*betak*rho;
            H = obj.F_all.*(1-Mk).*Delta;
            J = H'*H;
            constraints = abs(Delta) <= deltaN;

            evalc("optimize(constraints,J);");
            rho = value(rho);
            Delta = value(Delta);
        end
    end
end

