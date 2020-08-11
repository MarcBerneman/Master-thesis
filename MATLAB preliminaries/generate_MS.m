function [u,k_exc] = generate_MS(N,Fs,RMS,fmin,fmax,df,phases)
    f = (fmin:df:fmax).';
    fres = Fs/N;
    k_exc = round(f/fres);
    if nargin < 7
        phases = rand(size(k_exc))*2*pi;
    end
    U = zeros(N,1);
    U(k_exc+1) = exp(1j*phases);
    u = real(ifft(U));
    u = u*RMS/rms(u);
end

