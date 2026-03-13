function [TA, DYDA] = TAUA(fitParams, nDataPoints, MA, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth, measWavenumber)
% TAUA - Compute modeled transmittance and analytic derivatives
% Exact translation of FORTRAN SUBROUTINE TAUA
%
% FORTRAN source:
%   DO I=1,MA-2: IF(fitParams(I)<=0) fitParams(I)=1E-6  (enforce positivity)
%   DO K=1,nDataPoints: CALL CONV(measWavenumber(K),T1,fitParams)
%     TA(K)    = T1(1)
%     DYDA(K,I)= T1(I+1) for I=1..MA
%
% Inputs:
%   fitParams   - coefficient vector [p1..pNM, a1, b1] (MA x 1)
%   nDataPoints  - number of measured data points
%   MA   - total number of parameters (nMolecules+2)
%   optDepth   - optical depth array (nGridPoints+1 x nMolecules)
%   nGridPoints  - maximum grid index
%   nMolecules   - number of molecules
%   waveMin - minimum wavenumber (cm^-1)
%   fineGridStep  - wavenumber grid step (cm^-1)
%   ilsHalfWidth  - half spectral resolution (cm^-1)
%   measWavenumber   - measured wavenumber array (nDataPoints x 1)
%
% Outputs:
%   TA   - modeled transmittance (nDataPoints x 1)
%   DYDA - Jacobian matrix (nDataPoints x MA)
%          DYDA(K,I) = dTA(K)/dAY(I)

    % Initialize outputs
    TA   = zeros(nDataPoints, 1);
    DYDA = zeros(nDataPoints, MA);

    % FORTRAN: enforce positivity for molecule partial pressures
    % fitParams(I) <= 0 -> fitParams(I) = 1e-6  for I=1..nMolecules  (i.e., MA-2 molecules)
    for I = 1:MA-2
        if fitParams(I) <= 0.0
            fitParams(I) = 1.0e-6;
        end
    end

    % Loop over all measured data points
    for K = 1:nDataPoints
        V2 = measWavenumber(K);
        % CONV returns S(nMolecules+3): S(1)=T, S(2..nMolecules+1)=dT/dp_k, S(nMolecules+2)=dT/da1, S(nMolecules+3)=dT/db1
        T1 = CONV(V2, fitParams, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);
        TA(K) = T1(1);
        for I = 1:MA
            DYDA(K, I) = T1(I + 1);
        end
    end

end