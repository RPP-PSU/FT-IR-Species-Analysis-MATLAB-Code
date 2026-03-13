function OUTPT(WavenumberData, MeasuredTransmittance, CoefficientArray, ...
               CovarianceMatrix, ChiSquared, OutputFilePath1, OutputFilePath2, ...
               OutputFilePath3, NumMolecules, optDepth, nGridPoints, waveMin, fineGridStep, ilsHalfWidth, ...
               paramNames, nIterations, measWavenumber, measTransmittance, nDataPoints, nParams)
% OUTPT - Write output files and plot results
% Exact translation of FORTRAN SUBROUTINE OUTPT
%
% FORTRAN output structure:
%   Unit 21 (transmOutFile): measured vs calculated transmittance
%     FORMAT 101: (1X,F9.4,2(1X,F10.7),1X,E12.6)
%     Columns: measWavenumber(K), measTransmittance(K), TA(K), ERRT=measTransmittance(K)-TA(K)
%
%   Unit 22 (coeffOutFile): coefficients and covariance
%     FORMAT 102: (' CHI^2  =  ',E14.7)
%     FORMAT 103: (1X,A5,'  =  ',E12.6)      <- coefficient values
%     FORMAT 104: (1X,'SIG(',A5,') =  ',E12.6) <- std deviations
%     FORMAT 105: (5(1X,E12.6))               <- covariance matrix rows
%     FORMAT 110: (1X,'NO OF ITERATION IS    =',I5)
%
% FORTRAN FORMAT statements (exact):
%   101 FORMAT(1X,F9.4,2(1X,F10.7),1X,E12.6)
%   102 FORMAT(' CHI^2  =  ',E14.7)
%   103 FORMAT(1X,A5,'  =  ',E12.6)
%   104 FORMAT(1X,'SIG(',A5,') =  ',E12.6)
%   105 FORMAT(5(1X,E12.6))
%   110 FORMAT(1X,'NO OF ITERATION IS    =',I5)
%
% NOTE on FORTRAN E format:
%   FORTRAN E14.7 produces: -0.1234567E+02  (leading digit, 7 decimal places)
%   FORTRAN E12.6 produces: -0.123456E+02   (leading digit, 6 decimal places)
%   MATLAB fprintf %E uses: -1.234567E+02   (no leading zero convention)
%   Replication below uses sprintf with manual formatting to match FORTRAN
%   exponent style: 0.xxxxxxxE+xx
%
% Inputs:
%   WavenumberData       - measured wavenumber array (nDataPoints x 1)
%   MeasuredTransmittance- measured transmittance array (nDataPoints x 1)
%   CoefficientArray     - fitted coefficients fitParams (nParams x 1)
%   CovarianceMatrix     - covariance matrix covariance (nParams x nParams)
%   ChiSquared           - final chi-squared value
%   OutputFilePath1      - path for coeffOutFile (coefficients file)
%   OutputFilePath2      - path for transmOutFile (transmittance file)
%   OutputFilePath3      - unused (kept for interface compatibility)
%   NumMolecules         - nMolecules number of molecules
%   optDepth                   - optical depth array (nGridPoints+1 x nMolecules)
%   nGridPoints                  - maximum grid index
%   waveMin                 - minimum wavenumber (cm^-1)
%   fineGridStep                  - wavenumber grid step (cm^-1)
%   ilsHalfWidth                  - half spectral resolution (cm^-1)
%   paramNames                  - molecule name strings cell array (nParams x 1)
%   nIterations                  - number of iterations used
%   measWavenumber                   - wavenumber array from INDAT (nDataPoints x 1)
%   measTransmittance                   - measured transmittance from INDAT (nDataPoints x 1)
%   nDataPoints                  - number of data points
%   nParams                  - total number of parameters (nMolecules+2)

    % --- COMPUTE MODELED TRANSMITTANCE ---
    % Recompute TA using TAUA at final coefficients (same as FORTRAN OUTPT
    % which uses the TA array passed from the calling loop)
    fitParams = CoefficientArray;
    [TA, ~] = TAUA(fitParams, nDataPoints, nParams, optDepth, nGridPoints, NumMolecules, waveMin, fineGridStep, ilsHalfWidth, measWavenumber);

    % =========================================================
    % UNIT 21 -> transmOutFile: Transmittance comparison file
    % FORTRAN:
    %   DO K=1,nDataPoints
    %     ERRT = measTransmittance(K) - TA(K)
    %     WRITE(21,101) measWavenumber(K),measTransmittance(K),TA(K),ERRT
    %   101 FORMAT(1X,F9.4,2(1X,F10.7),1X,E12.6)
    % =========================================================
    fid21 = fopen(OutputFilePath2, 'w');
    if fid21 == -1
        error('OUTPT: cannot open output file: %s', OutputFilePath2);
    end

    for K = 1:nDataPoints
        ERRT = measTransmittance(K) - TA(K);
        % FORTRAN FORMAT 101: (1X,F9.4,2(1X,F10.7),1X,E12.6)
        % 1X  = 1 leading space
        % F9.4  = field width 9, 4 decimal places
        % 1X  = 1 space
        % F10.7 = field width 10, 7 decimal places  (x2)
        % 1X  = 1 space
        % E12.6 = field width 12, 6 decimal places scientific
        fprintf(fid21, ' %9.4f %10.7f %10.7f %s\n', ...
            measWavenumber(K), measTransmittance(K), TA(K), fortran_E12_6(ERRT));
    end
    fclose(fid21);

    % =========================================================
    % UNIT 22 -> coeffOutFile: Coefficients and covariance file
    % FORTRAN:
    %   WRITE(22,102) chiSq
    %   DO I=1,nParams: WRITE(22,103) paramNames(I),fitParams(I)
    %   DO I=1,nParams: WRITE(22,104) paramNames(I),covariance(I,I)
    %   DO I=1,nParams: WRITE(22,105) (covariance(I,J),J=1,nParams)
    %   WRITE(22,110) nIterations
    % =========================================================
    fid22 = fopen(OutputFilePath1, 'w');
    if fid22 == -1
        error('OUTPT: cannot open output file: %s', OutputFilePath1);
    end

    % FORMAT 102: ' CHI^2  =  ' followed by E14.7
    fprintf(fid22, ' CHI^2  =  %s\n', fortran_E14_7(ChiSquared));

    % FORMAT 103: (1X,A5,'  =  ',E12.6)
    % 1X = 1 space, A5 = 5-char name left-padded, '  =  ', E12.6
    for I = 1:nParams
        name_str = fortran_A5(paramNames{I});
        fprintf(fid22, ' %s  =  %s\n', name_str, fortran_E12_6(fitParams(I)));
    end

    % FORMAT 104: (1X,'SIG(',A5,') =  ',E12.6)
    for I = 1:nParams
        name_str = fortran_A5(paramNames{I});
        fprintf(fid22, ' SIG(%s) =  %s\n', name_str, fortran_E12_6(sqrt(abs(CovarianceMatrix(I,I)))));
    end

    % FORMAT 105: (5(1X,E12.6))
    % Each row of covariance matrix, 5 values per line with 1X separator
    for I = 1:nParams
        line_str = '';
        for J = 1:nParams
            line_str = [line_str, ' ', fortran_E12_6(CovarianceMatrix(I,J))]; %#ok<AGROW>
        end
        fprintf(fid22, '%s\n', line_str);
    end

    % FORMAT 110: (1X,'NO OF ITERATION IS    =',I5)
    fprintf(fid22, ' NO OF ITERATION IS    =%5d\n', nIterations);

    fclose(fid22);

    % =========================================================
    % PLOT: Measured vs Model Transmittance
    % Match existing plot format from current code exactly:
    %   Title:  'Measured vs. Model Transmittance'
    %           'ChiSquared: [value]'
    %   X-axis: 'Wavenumber (cm^{-1})'
    %   Y-axis: 'Transmittance'
    %   Blue solid:  Measured
    %   Red dashed:  Model
    %   Legend, grid on
    % =========================================================
    figure('Name', 'Measured vs. Model Transmittance', 'NumberTitle', 'off');
    plot(measWavenumber, measTransmittance, 'b-',  'DisplayName', 'Measured', 'LineWidth', 1.0);
    hold on;
    plot(measWavenumber, TA, 'r--', 'DisplayName', 'Model',    'LineWidth', 1.0);
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Transmittance');
    legend('show');
    title(sprintf('Measured vs. Model Transmittance\nChiSquared: %.7g', ChiSquared));
    grid on;
    hold off;

end


% =========================================================
% LOCAL HELPER FUNCTIONS
% Replicate FORTRAN E and F format output exactly
% =========================================================

function s = fortran_E14_7(val)
% Replicate FORTRAN E14.7 format
% FORTRAN E14.7: total width 14, 7 decimal places
% Form: [sign]0.dddddddE+ee  (leading zero before decimal)
% Example:  0.1234567E+02  or -0.1234567E-03
%
% FORTRAN sign convention:
%   positive: space before 0 (total 14 chars)
%   negative: minus before 0 (total 14 chars)

    if val == 0.0
        s = '  0.0000000E+00';
        s = s(end-13:end);   % trim to 14 chars
        return;
    end

    sgn = '';
    if val < 0
        sgn = '-';
        val = abs(val);
    else
        sgn = ' ';
    end

    exp10  = floor(log10(val));
    mant   = val / 10^exp10;       % 1.0 <= mant < 10
    mant   = mant / 10.0;          % 0.1 <= mant < 1.0  (FORTRAN leading zero)
    exp10  = exp10 + 1;

    % Round mantissa to 7 decimal places
    mant = round(mant * 1e7) / 1e7;

    % Handle carry from rounding
    if mant >= 1.0
        mant  = mant / 10.0;
        exp10 = exp10 + 1;
    end

    if exp10 >= 0
        s = sprintf('%s0.%07dE+%02d', sgn, round(mant * 1e7), exp10);
    else
        s = sprintf('%s0.%07dE-%02d', sgn, round(mant * 1e7), abs(exp10));
    end

    % Ensure total width = 14
    while length(s) < 14
        s = [' ', s]; %#ok<AGROW>
    end
end


function s = fortran_E12_6(val)
% Replicate FORTRAN E12.6 format
% FORTRAN E12.6: total width 12, 6 decimal places
% Form: [sign]0.ddddddE+ee
% Example:  0.123456E+02  or -0.123456E-03

    if val == 0.0
        s = '  0.000000E+00';
        s = s(end-11:end);   % trim to 12 chars
        return;
    end

    sgn = '';
    if val < 0
        sgn = '-';
        val = abs(val);
    else
        sgn = ' ';
    end

    exp10  = floor(log10(val));
    mant   = val / 10^exp10;
    mant   = mant / 10.0;
    exp10  = exp10 + 1;

    % Round mantissa to 6 decimal places
    mant = round(mant * 1e6) / 1e6;

    % Handle carry
    if mant >= 1.0
        mant  = mant / 10.0;
        exp10 = exp10 + 1;
    end

    if exp10 >= 0
        s = sprintf('%s0.%06dE+%02d', sgn, round(mant * 1e6), exp10);
    else
        s = sprintf('%s0.%06dE-%02d', sgn, round(mant * 1e6), abs(exp10));
    end

    % Ensure total width = 12
    while length(s) < 12
        s = [' ', s]; %#ok<AGROW>
    end
end


function s = fortran_A5(name)
% Replicate FORTRAN A5 format
% A5 = exactly 5 characters, left-justified, space-padded on right
% FORTRAN paramNames array is CHARACTER*8 but A5 truncates/pads to 5 chars

    name = strtrim(name);
    if length(name) > 5
        name = name(1:5);
    end
    % Left-justify with trailing spaces to width 5
    s = sprintf('%-5s', name);
end