function V = RLFR(NV1, waveMin, fineGridStep)
% RLFR - Convert array index to wavenumber
% Exact translation of FORTRAN FUNCTION RLFR
%
% FORTRAN source:
%   FUNCTION RLFR(NV1)
%   RLFR = FLOAT(NV1 - 1)*fineGridStep + waveMin - 30.D0
%
% Inputs:
%   NV1  - array index (1-based integer)
%   waveMin - minimum wavenumber of spectral range (cm^-1)
%   fineGridStep  - wavenumber grid step (cm^-1)
%
% Output:
%   V    - wavenumber (cm^-1) corresponding to index NV1

    V = double(NV1 - 1) * fineGridStep + waveMin - 30.0;
end
