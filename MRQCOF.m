function [hessian, gradient, chiSq] = MRQCOF(TA_dummy, fitParams, nDataPoints, MA, measWavenumber, measTransmittance, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth)
% MRQCOF - Compute curvature matrix hessian and vector gradient for LM method
% Exact translation of FORTRAN SUBROUTINE MRQCOF
%
% FORTRAN source:
%   All parameters IA(J)=1 (all parameters are free - no fixed params)
%   MFIT = MA (all parameters fitted)
%   hessian(j,k) accumulated from DYDA products
%   gradient(j)    accumulated from residual * DYDA
%   chiSq      sum of squared residuals
%   Symmetry: hessian(k,j) = hessian(j,k) filled at end
%
% Inputs:
%   TA_dummy - unused (matches FORTRAN call signature)
%   fitParams       - current coefficient vector (MA x 1)
%   nDataPoints      - number of measured data points
%   MA       - total parameters (nMolecules+2)
%   measWavenumber       - measured wavenumber array (nDataPoints x 1)
%   measTransmittance       - measured transmittance array (nDataPoints x 1)
%   optDepth       - optical depth array (nGridPoints+1 x nMolecules)
%   nGridPoints      - maximum grid index
%   nMolecules       - number of molecules
%   waveMin     - minimum wavenumber (cm^-1)
%   fineGridStep      - wavenumber grid step (cm^-1)
%   ilsHalfWidth      - half spectral resolution (cm^-1)
%
% Outputs:
%   hessian - curvature matrix (MA x MA)
%   gradient  - derivative vector (MA x 1)
%   chiSq - chi-squared value (sum of squared residuals)

    % FORTRAN: all IA(J)=1, so MFIT=MA
    MFIT = MA;

    % Initialize
    hessian = zeros(MA, MA);
    gradient  = zeros(MA, 1);
    chiSq = 0.0;

    % Compute modeled transmittance and Jacobian via TAUA
    [TA, DYDA] = TAUA(fitParams, nDataPoints, MA, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth, measWavenumber);

    % Accumulate hessian and gradient
    % FORTRAN loop structure (lower triangle only, then symmetrize):
    %   DO I=1,nDataPoints
    %     DY = measTransmittance(I) - TA(I)
    %     DO L=1,MA (j index, only if IA(L)!=0)
    %       DO M=1,L (k index, only if IA(M)!=0)
    %         hessian(j,k) += DYDA(I,L)*DYDA(I,M)
    %       gradient(j) += DY*DYDA(I,L)
    %     chiSq += DY*DY

    for I = 1:nDataPoints
        DY = measTransmittance(I) - TA(I);
        for L = 1:MFIT
            for M = 1:L
                hessian(L, M) = hessian(L, M) + DYDA(I, L) * DYDA(I, M);
            end
            gradient(L) = gradient(L) + DY * DYDA(I, L);
        end
        chiSq = chiSq + DY * DY;
    end

    % Symmetrize: hessian(k,j) = hessian(j,k) for j=2..MFIT, k=1..j-1
    for J = 2:MFIT
        for K = 1:J-1
            hessian(K, J) = hessian(J, K);
        end
    end

end