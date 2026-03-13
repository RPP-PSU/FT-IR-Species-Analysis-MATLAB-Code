function A = APND(X1, ilsHalfWidth)
% APND - Instrument Line Shape (sinc-squared) function
% Exact translation of FORTRAN FUNCTION APND
%
% FORTRAN source:
%   DATA PI1,XNL/ 0.314159D+1, 2.479D+19 /  <- PI1 = 3.14159 exactly
%   ilsHalfSpan = 0.5D0/ilsHalfWidth
%   IF(X1.EQ.0.) A1=1 ELSE A1=SIN(PI*X1*ilsHalfSpan)/(X1*PI*ilsHalfSpan)
%   APND = A1*A1*ilsHalfSpan
%
% Inputs:
%   X1  - wavenumber offset (scalar)
%   ilsHalfWidth - half spectral resolution of FTIR instrument (cm^-1)
%
% Output:
%   A   - ILS value at offset X1

    PI1  = 0.314159D+1;    % FORTRAN exact constant: 3.14159
    ilsHalfSpan = 0.5 / ilsHalfWidth;

    if X1 == 0.0
        A1 = 1.0;
    else
        A1 = sin(PI1 * X1 * ilsHalfSpan) / (X1 * PI1 * ilsHalfSpan);
    end

    A = A1 * A1 * ilsHalfSpan;
end
