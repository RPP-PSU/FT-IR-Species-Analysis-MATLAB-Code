function NV = INFR(VI, waveMin, fineGridStep)
% INFR - Convert wavenumber to array index
% Exact translation of FORTRAN FUNCTION INFR
%
% FORTRAN source:
%   FUNCTION INFR(VI)
%   INFR = DINT((VI - waveMin + 30.D0)/fineGridStep) + 1
%
% DINT in FORTRAN truncates toward zero (equivalent to fix() in MATLAB)
% NOT round() - this distinction is critical for index correctness
%
% Inputs:
%   VI   - wavenumber value (cm^-1)
%   waveMin - minimum wavenumber of spectral range (cm^-1)
%   fineGridStep  - wavenumber grid step (cm^-1)
%
% Output:
%   NV   - integer array index (1-based)

    NV = fix((VI - waveMin + 30.0) / fineGridStep) + 1;
end
