function optDepth = OPTD(V, lineStrength, halfWidth, KK, optDepth, nGridPoints, nMolecules, waveMin, waveMax, fineGridStep, temperature, pathLength, IFLG)
% OPTD - Compute optical depth contribution from a single spectral line
% Exact translation of FORTRAN SUBROUTINE OPTD
%
% FORTRAN source:
%   DATA PI1,XNL/ 0.314159D+1, 2.479D+19 /
%   XN = XNL*296.D0/temperature
%   WVLW = waveMin - 15.D0
%   WVHG = waveMax + 15.D0
%   nConvPoints  = DINT(15.D0/fineGridStep)
%
%   IF IFLG==1: zero optDepth(1..nGridPoints+1, 1..nMolecules), set IFLG=2
%
%   NV = INFR(V)   <- DINT-based index
%   Three cases for NTM, IMAX, VI based on V vs WVLW/WVHG
%
%   XKI = pathLength*XN/PI1*lineStrength*halfWidth
%   DO J=1,IMAX: optDepth(NTM,K) += XKI/((VI-V)^2 + halfWidth^2)
%
% CRITICAL: IFLG logic
%   IFLG=1 on entry to first OPTD call per scan -> zero optDepth, set IFLG=2
%   IFLG=2 on subsequent calls -> only accumulate
%
% Inputs:
%   V    - line center wavenumber (cm^-1)
%   lineStrength  - temperature-corrected line strength
%   halfWidth - half-width at half-maximum (cm^-1)
%   KK   - molecule index (1-based column of optDepth)
%   optDepth   - optical depth array (nGridPoints+1 x nMolecules), modified in place
%   nGridPoints  - maximum grid index = fix((waveMax-waveMin+60)/fineGridStep)
%   nMolecules   - number of molecules
%   waveMin - minimum wavenumber of spectral range (cm^-1)
%   waveMax - maximum wavenumber of spectral range (cm^-1)
%   fineGridStep  - wavenumber grid step (cm^-1)
%   temperature   - temperature (K)
%   pathLength   - path length (cm)
%   IFLG - initialization flag (1=first call, 2=subsequent calls)
%          NOTE: caller must check returned IFLG and pass it back in
%
% Outputs:
%   optDepth   - updated optical depth array

    % FORTRAN DATA statement constants - exact values
    PI1 = 0.314159D+1;    % 3.14159 exactly as in FORTRAN
    XNL = 2.479D+19;      % Loschmidt number approximation

    % Number density at temperature temperature (FORTRAN: XN = XNL*296/temperature)
    XN = XNL * 296.0 / temperature;

    % Wing boundaries (FORTRAN: WVLW=waveMin-15, WVHG=waveMax+15)
    WVLW = waveMin - 15.0;
    WVHG = waveMax + 15.0;

    % FORTRAN: nConvPoints = DINT(15.D0/fineGridStep)  <- truncation via fix()
    nConvPoints = fix(15.0 / fineGridStep);

    % --- IFLG INITIALIZATION ---
    % FORTRAN: IF(IFLG.EQ.1) zero optDepth, set IFLG=2
    % This happens on the FIRST call to OPTD for each scan
    if IFLG == 1
        for J = 1:nGridPoints+1
            for K = 1:nMolecules
                optDepth(J, K) = 0.0;
            end
        end
        IFLG = 2;  % caller must use returned optDepth; IFLG state tracked in INPUT.m
    end

    % Molecule column index
    K = KK;

    % --- FIND CENTER INDEX ---
    % FORTRAN: NV = INFR(V) using DINT
    NV = INFR(V, waveMin, fineGridStep);

    % --- DETERMINE WINDOW (three cases from FORTRAN) ---
    if V < WVLW
        % Line center below lower wing boundary
        NTM  = 1;
        IMAX = NV + nConvPoints;
        VI   = waveMin - 30.0;    % FORTRAN: VI = waveMin - 30.D0
    elseif V > WVHG
        % Line center above upper wing boundary
        NTM  = NV - nConvPoints;
        IMAX = nGridPoints + nConvPoints + 2 - NV;
        VI   = RLFR(NTM, waveMin, fineGridStep);
    else
        % Line center within [WVLW, WVHG] (standard case)
        NTM  = NV - nConvPoints;
        IMAX = 2 * nConvPoints + 1;
        VI   = RLFR(NTM, waveMin, fineGridStep);
    end

    % --- LORENTZIAN PREFACTOR ---
    % FORTRAN: XKI = pathLength*XN/PI1*lineStrength*halfWidth
    XKI   = pathLength * XN / PI1 * lineStrength * halfWidth;
    HWHM2 = halfWidth * halfWidth;

    % --- ACCUMULATE OPTICAL DEPTH ---
    % FORTRAN: DO J=1,IMAX
    %   XK = XKI/((VI-V)^2 + HWHM2)
    %   optDepth(NTM,K) += XK
    %   VI += fineGridStep; NTM += 1
    for J = 1:IMAX
        if NTM > nGridPoints + 1
            break;    % safety: do not exceed array bounds
        end
        if NTM >= 1
            XK          = XKI / ((VI - V)^2 + HWHM2);
            optDepth(NTM, K)  = optDepth(NTM, K) + XK;
        end
        VI  = VI  + fineGridStep;
        NTM = NTM + 1;
    end

end
