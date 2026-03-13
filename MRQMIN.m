function [fitParams, covariance, hessian, chiSq, lambdaDamp, stepImproved, OCHISQ_out, ATRY_out, BETA_out, DA_out] = ...
    MRQMIN(TA, nDataPoints, fitParams, MA, covariance, hessian, chiSq, lambdaDamp, stepImproved, ...
           measWavenumber, measTransmittance, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth, ...
           OCHISQ_in, ATRY_in, BETA_in, DA_in)
% MRQMIN - Levenberg-Marquardt nonlinear least squares controller
% Exact translation of FORTRAN SUBROUTINE MRQMIN
%
% FORTRAN source logic:
%   1. If lambdaDamp < 0 (first call):
%      - Set lambdaDamp=0.001
%      - Call MRQCOF to initialize chiSq, hessian, gradient
%      - Save chiSqPrev=chiSq, paramsTrial=fitParams
%   2. Build trial covariance from hessian with diagonal augmented by (1+lambdaDamp)
%   3. Call GAUSSJ to solve for paramStep
%   4. If lambdaDamp==0: return (final covariance call)
%   5. Try paramsTrial = fitParams + paramStep  (with positivity clamp for molecule params)
%   6. Call MRQCOF with paramsTrial to get new chiSq
%   7. If improved: stepImproved=0, lambdaDamp*=0.1, update hessian/gradient/fitParams
%      Else:        stepImproved=1, lambdaDamp*=10,  restore chiSq=chiSqPrev
%
% Saved state (replaces FORTRAN SAVE statement):
%   chiSqPrev, paramsTrial, gradient, paramStep passed in and out explicitly
%
% stepImproved=0 means solution improved this step
% stepImproved=1 means solution did not improve this step
%
% POSITIVITY CLAMP:
%   After building paramsTrial = fitParams + paramStep, indices 1..MA-2 (molecule partial
%   pressures) are clamped to >= 1e-10.
%   Negative partial pressures are physically impossible and cause
%   exp(-optDepth*fitParams) -> exp(+large) which makes the model diverge.
%   This clamp replicates the positivity enforcement in TAUA.m but
%   applied to the TRIAL vector before MRQCOF evaluates it.

    %#ok<*INUSL>

    MFIT = MA;

    % =========================================================
    % PRE-ASSIGN ALL OUTPUTS TO SAFE DEFAULTS
    % Guarantees no "output not assigned" error on any code path
    % =========================================================
    OCHISQ_out = OCHISQ_in;
    ATRY_out   = ATRY_in;
    BETA_out   = BETA_in;
    DA_out     = DA_in;

    % =========================================================
    % INITIALIZATION CALL  (lambdaDamp < 0)
    % FORTRAN: IF(lambdaDamp.LT.0) THEN
    % =========================================================
    if lambdaDamp < 0.0

        % FORTRAN: lambdaDamp = 0.001
        lambdaDamp = 0.001;

        % FORTRAN: CALL MRQCOF(fitParams,hessian,gradient,chiSq)
        [hessian, gradient, chiSq] = MRQCOF([], fitParams, nDataPoints, MA, measWavenumber, measTransmittance, ...
            optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);

        % FORTRAN: chiSqPrev = chiSq
        OCHISQ_out = chiSq;

        % FORTRAN: DO J=1,MA: paramsTrial(J)=fitParams(J)
        ATRY_out = fitParams;
        BETA_out = gradient;
        DA_out   = zeros(MA, 1);

        % Return current hessian as covariance on init
        covariance = hessian;
        stepImproved    = 0;
        return;
    end

    % Recover saved state from previous call
    chiSqPrev = OCHISQ_in;
    paramsTrial   = ATRY_in;
    gradient   = BETA_in;
    paramStep     = DA_in;

    % =========================================================
    % BUILD TRIAL COVARIANCE MATRIX (augmented diagonal)
    % FORTRAN:
    %   DO j=1,mfit
    %     DO k=1,mfit: covariance(j,k) = hessian(j,k)
    %     covariance(j,j) = hessian(j,j)*(1+lambdaDamp)
    %     paramStep(j)      = gradient(j)
    % =========================================================
    j = 0;
    for l = 1:MA
        j = j + 1;
        k = 0;
        for m = 1:MA
            k = k + 1;
            covariance(j, k) = hessian(j, k);
        end
        covariance(j, j) = hessian(j, j) * (1.0 + lambdaDamp);
        paramStep(j)       = gradient(j);
    end

    % =========================================================
    % SOLVE LINEAR SYSTEM VIA GAUSSJ
    % FORTRAN: CALL GAUSSJ(covariance, mfit, paramStep, 1)
    % =========================================================
    DA_mat = paramStep(:);
    [covariance, DA_mat] = GAUSSJ(covariance, MFIT, DA_mat, 1);
    paramStep = DA_mat(:);

    % =========================================================
    % FINAL COVARIANCE CALL  (lambdaDamp == 0)
    % FORTRAN: IF(lambdaDamp.EQ.0) RETURN
    % All outputs pre-assigned above - safe to return
    % =========================================================
    if lambdaDamp == 0.0
        OCHISQ_out = chiSqPrev;
        ATRY_out   = paramsTrial;
        BETA_out   = gradient;
        DA_out     = paramStep;
        stepImproved         = 0;
        return;
    end

    % =========================================================
    % BUILD TRIAL PARAMETER VECTOR
    % FORTRAN: DO l=1,ma: paramsTrial(l) = fitParams(l) + paramStep(j)
    % =========================================================
    j = 0;
    for l = 1:MA
        j = j + 1;
        paramsTrial(l) = fitParams(l) + paramStep(j);
    end

    % =========================================================
    % POSITIVITY CLAMP FOR MOLECULE PARTIAL PRESSURES
    % Indices 1..MA-2 are molecule concentrations (MA-2 = nMolecules)
    % Indices MA-1 and MA are baseline a1 and b1 - not clamped
    %
    % Negative partial pressures are physically impossible.
    % Without this clamp, a negative paramsTrial(k) causes:
    %   exp(-optDepth * paramsTrial(k)) -> exp(+large) -> model transmittance >> 1
    % which makes MRQCOF return a huge chiSq and the LM fitter
    % diverges or gets stuck in an unphysical region.
    %
    % Clamp to 1e-10 (not zero) to keep the Jacobian well-defined.
    % =========================================================
    for l = 1:MA-2
        if paramsTrial(l) <= 0.0
            paramsTrial(l) = 1.0e-10;
        end
    end

    % =========================================================
    % EVALUATE TRIAL SOLUTION
    % FORTRAN: CALL MRQCOF(paramsTrial,covariance,paramStep,CHISQ_trial)
    % =========================================================
    [COVAR_trial, DA_trial, CHISQ_trial] = MRQCOF([], paramsTrial, nDataPoints, MA, ...
        measWavenumber, measTransmittance, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);

    % =========================================================
    % ACCEPT OR REJECT TRIAL SOLUTION
    % =========================================================
    if CHISQ_trial < chiSqPrev

        % --- SOLUTION IMPROVED ---
        % FORTRAN: stepImproved=0, lambdaDamp*=0.1, update hessian/gradient/fitParams
        stepImproved     = 0;
        lambdaDamp = 0.1 * lambdaDamp;
        chiSqPrev = CHISQ_trial;
        chiSq  = CHISQ_trial;

        % FORTRAN: hessian = COVAR_trial, gradient = DA_trial, fitParams = paramsTrial
        j = 0;
        for l = 1:MA
            j = j + 1;
            k = 0;
            for m = 1:MA
                k = k + 1;
                hessian(j, k) = COVAR_trial(j, k);
            end
            gradient(j) = DA_trial(j);
            fitParams(l)   = paramsTrial(l);
        end
        covariance = COVAR_trial;

    else

        % --- SOLUTION DID NOT IMPROVE ---
        % FORTRAN: stepImproved=1, lambdaDamp*=10, restore chiSq=chiSqPrev
        stepImproved     = 1;
        lambdaDamp = 10.0 * lambdaDamp;
        chiSq  = chiSqPrev;

    end

    % Pass saved state forward to next call
    OCHISQ_out = chiSqPrev;
    ATRY_out   = paramsTrial;
    BETA_out   = gradient;
    DA_out     = paramStep;

end