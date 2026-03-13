function S = CONV(V1, fitParams, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth)
% CONV - Convolution of true transmittance with instrument line shape
% Exact translation of FORTRAN SUBROUTINE CONV
%
% TWO SEPARATE GRID PARAMETERS:
%   fineGridStep - fine wavenumber grid step (cm^-1) used to build optDepth array
%         FORTRAN input line 1: 5D0 1.D0 fineGridStep
%         Typical value: 0.015 cm^-1
%         Controls: nConvPoints = fix(30/fineGridStep), grid index arithmetic
%
%   ilsHalfWidth - half-width of sinc^2 ILS (cm^-1), matches data point spacing
%         FORTRAN input line 2: waveMin waveMax ilsHalfWidth
%         Typical value: 1.28565 cm^-1
%         Controls: APND(X1, ilsHalfWidth) -> ilsHalfSpan = 0.5/ilsHalfWidth
%
% FORTRAN trapezoidal integration logic:
%   - Endpoint half-weights applied ONLY when endpoint is OUT OF BOUNDS
%   - Interior points (gridIdxLeft to gridIdxRight inclusive) always use full weight
%   - Final result scaled by fineGridStep (grid step)
%
% Inputs:
%   V1   - measurement wavenumber (cm^-1)
%   fitParams   - coefficient vector [p1..pNM, a1, b1] (nMolecules+2 x 1)
%   optDepth   - optical depth array (nGridPoints+1 x nMolecules), built on fine fineGridStep grid
%   nGridPoints  - maximum fine grid index = fix((waveMax-waveMin+60)/fineGridStep)
%   nMolecules   - number of molecules
%   waveMin - minimum wavenumber (cm^-1)
%   fineGridStep  - fine optDepth grid step (cm^-1) — NOT the data spacing
%   ilsHalfWidth  - ILS half-width (cm^-1) — matches data point spacing
%
% Output:
%   S    - result vector (nMolecules+3 x 1):
%          S(1)       = modeled transmittance T(V1)
%          S(2..nMolecules+1) = dT/d(fitParams(k)) for k=1..nMolecules
%          S(nMolecules+2)    = dT/d(a1)
%          S(nMolecules+3)    = dT/d(b1)

    % FORTRAN: nConvPoints = DINT(30.D0/fineGridStep)
    % Uses FINE grid fineGridStep — gives ~2000 steps with fineGridStep=0.015
    nConvPoints = fix(30.0 / fineGridStep);

    NA = nMolecules + 1;   % baseline a1 index in fitParams
    NB = nMolecules + 2;   % baseline b1 index in fitParams

    % Find center grid index using INFR (DINT-based, matches FORTRAN)
    NV  = INFR(V1, waveMin, fineGridStep);
    gridIdxLeft = NV - nConvPoints;
    gridIdxRight = NV + nConvPoints;

    % Initialize output S vector
    S = zeros(nMolecules + 3, 1);

    % =========================================================
    % LEFT ENDPOINT
    % Half-weight applied ONLY if gridIdxLeft is out of bounds (< 1)
    % If gridIdxLeft >= 1 it is valid and included in full-weight interior loop
    % =========================================================
    if gridIdxLeft < 1
        gridIdxLeft = 1;
    else
        V3   = RLFR(gridIdxLeft, waveMin, fineGridStep);
        totalOptDepth = 0.0;
        for K = 1:nMolecules
            totalOptDepth = totalOptDepth - optDepth(gridIdxLeft, K) * fitParams(K);
        end
        % APND uses ilsHalfWidth (ILS half-width) — NOT fineGridStep
        A2 = exp(totalOptDepth) * APND(V1 - V3, ilsHalfWidth);

        S(1) = S(1) - 0.5 * (fitParams(NA) + fitParams(NB) * V3) * A2;
        for K = 1:nMolecules
            S(K+1) = S(K+1) + 0.5 * A2 * optDepth(gridIdxLeft,K) * (fitParams(NA) + fitParams(NB)*V3);
        end
        S(nMolecules+2) = S(nMolecules+2) - 0.5 * A2;
        S(nMolecules+3) = S(nMolecules+3) - 0.5 * A2 * V3;
    end

    % =========================================================
    % RIGHT ENDPOINT
    % Half-weight applied ONLY if gridIdxRight is out of bounds (> nGridPoints+1)
    % =========================================================
    if gridIdxRight > nGridPoints + 1
        gridIdxRight = nGridPoints + 1;
    else
        V3   = RLFR(gridIdxRight, waveMin, fineGridStep);
        totalOptDepth = 0.0;
        for K = 1:nMolecules
            totalOptDepth = totalOptDepth - optDepth(gridIdxRight, K) * fitParams(K);
        end
        % APND uses ilsHalfWidth (ILS half-width) — NOT fineGridStep
        A2 = exp(totalOptDepth) * APND(V1 - V3, ilsHalfWidth);

        S(1) = S(1) - 0.5 * (fitParams(NA) + fitParams(NB) * V3) * A2;
        for K = 1:nMolecules
            S(K+1) = S(K+1) + 0.5 * A2 * optDepth(gridIdxRight,K) * (fitParams(NA) + fitParams(NB)*V3);
        end
        S(nMolecules+2) = S(nMolecules+2) - 0.5 * A2;
        S(nMolecules+3) = S(nMolecules+3) - 0.5 * A2 * V3;
    end

    % =========================================================
    % INTERIOR LOOP: gridIdxLeft to gridIdxRight inclusive, full weight
    % FORTRAN: DO J = gridIdxLeft,gridIdxRight
    % =========================================================
    for J = gridIdxLeft:gridIdxRight
        V3   = RLFR(J, waveMin, fineGridStep);
        DV   = V1 - V3;
        totalOptDepth = 0.0;
        for K = 1:nMolecules
            totalOptDepth = totalOptDepth - optDepth(J, K) * fitParams(K);
        end
        % APND uses ilsHalfWidth (ILS half-width) — NOT fineGridStep
        A2 = exp(totalOptDepth) * APND(DV, ilsHalfWidth);

        S(1) = S(1) + (fitParams(NA) + fitParams(NB) * V3) * A2;
        for K = 1:nMolecules
            S(K+1) = S(K+1) - A2 * optDepth(J,K) * (fitParams(NA) + fitParams(NB)*V3);
        end
        S(nMolecules+2) = S(nMolecules+2) + A2;
        S(nMolecules+3) = S(nMolecules+3) + A2 * V3;
    end

    % =========================================================
    % SCALE BY fineGridStep — trapezoidal rule grid step multiplier
    % FORTRAN: DO K=1,nMolecules+3: S(K) = fineGridStep * S(K)
    % Uses FINE grid fineGridStep
    % =========================================================
    S = fineGridStep .* S;

end