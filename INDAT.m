function [measWavenumber, measTransmittance, NDA, nDataPoints] = INDAT(scanFilePath, waveMin, waveMax)
% INDAT - Read measured FTIR transmittance data from file
% Exact translation of FORTRAN SUBROUTINE INDAT
%
% FORTRAN source:
%   OPEN(4, FILE=scanFilePath, STATUS='OLD')
%   NDT  = 2451                        <- hardcoded total lines to read
%   NDA  = 0
%   nDataPoints  = 0
%   DO I=1,NDT
%     READ(4,*) VTP,TMTP              <- free-format read: 2 columns
%     IF(VTP>=waveMin .AND. VTP<=waveMax)  <- wavenumber range filter
%       nDataPoints = nDataPoints+1
%       measWavenumber(nDataPoints) = VTP
%       measTransmittance(nDataPoints) = TMTP
%   CLOSE(4)
%
%   <- REVERSAL BLOCK (FORTRAN: reorder to ascending wavenumber) ->
%   DO i=1,nDataPoints
%     dummy1(nDataPoints+1-i) = measWavenumber(i)
%     dummy2(nDataPoints+1-i) = measTransmittance(i)
%   DO i=1,nDataPoints
%     measWavenumber(i) = dummy1(i)
%     measTransmittance(i) = dummy2(i)
%
% KEY POINTS:
%   1. NDT=2451 is hardcoded exactly as in FORTRAN - not dynamic
%   2. Data files are in DECREASING wavenumber order (as confirmed by
%      NO_N2O.1 sample: 3999.64 down to 500.12)
%   3. FORTRAN reversal converts to ASCENDING order for fitting
%   4. NDA is initialized to 0 and never modified in FORTRAN INDAT
%      (NDA is set elsewhere via COMMON /TMEAS/ - kept here for
%       compatibility with COMMON block equivalents)
%   5. Free-format READ(*) handles any whitespace between columns
%      -> use textscan() to match FORTRAN behavior
%   6. measWavenumber and measTransmittance are dimensioned to 1000 in FORTRAN COMMON /TMEAS/
%      -> enforce max 1000 data points
%
% Inputs:
%   scanFilePath - full path to measured FTIR data file
%   waveMin   - minimum wavenumber of spectral range (cm^-1)
%   waveMax   - maximum wavenumber of spectral range (cm^-1)
%
% Outputs:
%   measWavenumber  - wavenumber array, ascending order (nDataPoints x 1)
%   measTransmittance  - transmittance array, ascending order (nDataPoints x 1)
%   NDA - initialized to 0 (matches FORTRAN COMMON /TMEAS/ NDA=0)
%   nDataPoints - number of data points within [waveMin, waveMax]

    % FORTRAN hardcoded constants
    NDT     = 2451;     % exact FORTRAN: NDT=2451
    MAXPTS  = 1000;     % FORTRAN DIMENSION measWavenumber(1000),measTransmittance(1000)

    % Initialize outputs matching FORTRAN
    measWavenumber  = zeros(MAXPTS, 1);
    measTransmittance  = zeros(MAXPTS, 1);
    NDA = 0;             % FORTRAN: NDA=0, set in COMMON /TMEAS/
    nDataPoints = 0;             % FORTRAN: nDataPoints=0

    % --- OPEN FILE ---
    % FORTRAN: OPEN(4,FILE=scanFilePath,STATUS='OLD')
    fid = fopen(scanFilePath, 'r');
    if fid == -1
        error('INDAT: cannot open measured data file: %s', scanFilePath);
    end

    % --- READ EXACTLY NDT=2451 LINES ---
    % FORTRAN: DO I=1,NDT: READ(4,*) VTP,TMTP
    % Free-format READ(*) handles any whitespace/tab delimiter
    % textscan replicates FORTRAN free-format read behavior
    raw = textscan(fid, '%f %f', NDT, 'CollectOutput', true);
    fclose(fid);

    if isempty(raw) || isempty(raw{1})
        error('INDAT: no data read from file: %s', scanFilePath);
    end

    raw_data = raw{1};    % NDT x 2 matrix: col1=VTP, col2=TMTP
    nread    = size(raw_data, 1);

    if nread < NDT
        warning('INDAT: file %s has %d lines, expected %d (NDT=2451)', ...
            scanFilePath, nread, NDT);
    end

    % --- APPLY WAVENUMBER RANGE FILTER ---
    % FORTRAN: IF(VTP>=waveMin .AND. VTP<=waveMax) THEN nDataPoints=nDataPoints+1
    for I = 1:nread
        VTP  = raw_data(I, 1);
        TMTP = raw_data(I, 2);
        if VTP >= waveMin && VTP <= waveMax
            nDataPoints        = nDataPoints + 1;
            measWavenumber(nDataPoints)    = VTP;
            measTransmittance(nDataPoints)    = TMTP;
        end
    end

    % Guard: nDataPoints must not exceed FORTRAN array dimension of 1000
    if nDataPoints > MAXPTS
        error('INDAT: nDataPoints=%d exceeds FORTRAN array limit of %d', nDataPoints, MAXPTS);
    end

    % Trim to nDataPoints
    measWavenumber = measWavenumber(1:nDataPoints);
    measTransmittance = measTransmittance(1:nDataPoints);

    % --- REVERSAL TO ASCENDING WAVENUMBER ORDER ---
    % FORTRAN:
    %   DO i=1,nDataPoints: dummy1(nDataPoints+1-i)=measWavenumber(i); dummy2(nDataPoints+1-i)=measTransmittance(i)
    %   DO i=1,nDataPoints: measWavenumber(i)=dummy1(i); measTransmittance(i)=dummy2(i)
    % This reverses the array in place using two explicit loops
    % Net effect: measWavenumber and measTransmittance go from descending to ascending order
    dummy1 = zeros(nDataPoints, 1);
    dummy2 = zeros(nDataPoints, 1);

    for i = 1:nDataPoints
        dummy1(nDataPoints + 1 - i) = measWavenumber(i);
        dummy2(nDataPoints + 1 - i) = measTransmittance(i);
    end

    for i = 1:nDataPoints
        measWavenumber(i) = dummy1(i);
        measTransmittance(i) = dummy2(i);
    end

end