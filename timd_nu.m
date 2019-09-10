function [TIMDr,TIMDf] = timd_nu(Y,f1,f2,fs,f,nmax,Nsl)

% This function computes the TIMD for a given input spectra.
% 
% Input Parameters
% Y      - Absolute value of the FFT response. Only provide the half of the
%           result from FFT
% f1     - Modulation frequency
% f2     - Fundamental frequency
% fs     - Sampling frequency
% f       - Frequency vector
% nmax - Maximum number of harmonics to consider
% Nsl    - Number of side lobes around the main lobe to consider. Specify
%           if the signal is windowed.
%
% Output Parameters
% TIMDr - TIMD excluding the fundamental in the denominator
% TIMDf - TIMD including the fundamental in the denominator

% compute the harmonics below f2
negf2 = 0;
for i = 1:nmax
    fh = f2-i*f1;
    % Break the loop if the harmonic frequency is less than 0
    if (fh < 0)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    negf2 = negf2 + sum(Y(ind(1:Nsl)).^2);
end

% compute the harmonics above f2
posf2 = 0;
for i = 1:nmax
    fh = f2+i*f1;
    % Break the loop if the harmonic frequency is greater than fs/2
    if (fh > fs/2)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    posf2 = posf2 + sum(Y(ind(1:Nsl)).^2);
end

% compute the value at f2
[~,ind] = sort(abs(f-f2));
atf2 = sum(Y(ind(1:Nsl)).^2);

% Now build the numerator for both TIMDr and TIMDf
num = negf2+posf2;

% Denominator for TIMDr
denr = negf2+atf2+posf2;

% Denominator for TIMDf
denf = atf2;

% Compute TIMDr
TIMDr = 100*sqrt(num/denr);

% Compute TIMDf
TIMDf = 100*sqrt(num/denf);
end