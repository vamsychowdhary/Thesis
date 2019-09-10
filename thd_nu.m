function [THDr,THDf] = thd_nu(Y,f0,f,nmax,Nsl)

% This function computes the THD for a given input spectra.
% 
% Input Parameters
% Y      - Absolute value of the FFT response. Only provide the half of the
%           result from FFT
% f0     - Fundamental frequency
% f       - Frequency vector
% nmax - Maximum number of harmonics to consider
% Nsl    - Number of side lobes around the main lobe to consider. Specify
%           if the signal is windowed.
%
% Output Parameters
% THDr - THD excluding the fundamental in the denominator
% THDf - THD including the fundamental in the denominator

% compute the numerator for THDr and THDf
num = 0;
for i = 2:nmax
    fh = i*f0;
    % Break the loop if order exceeds maximum frequency available
    if (fh > max(f))
        break;
    end
    [~,ind] = sort(abs(f-fh));
    num = num + sum(Y(ind(1:Nsl)).^2);
end

% compute the denominator for THDr
denr = 0;
for i = 1:nmax
    fh = i*f0;
    % Break the loop if order exceeds maximum frequency available    
    if (fh > max(f))
        break;
    end
    [~,ind] = sort(abs(f-fh));
    denr = denr + sum(Y(ind(1:Nsl)).^2);
end

% compute the denominator for THDf
[~,ind] = sort(abs(f-f0));
denf = sum(Y(ind(1:Nsl)).^2);

% Compute THDr
THDr = 100*sqrt(num/denr);

% Compute THDf
THDf = 100*sqrt(num/denf);
end