function QT = QofT_matlab(Tgrid, Qgrid, Tval)
% QofT_matlab - Partition function interpolation
% Exact translation of FORTRAN SUBROUTINE QofT
%
% FORTRAN source:
%   DATA Tref,eps / 296.d0, 1.D-05/
%   ivec = isovec(Mol) + iso       <- handled by caller, ivec row passed as Qgrid
%   Itemp = DINT(Tout) - 295       <- DINT = truncation toward zero = fix()
%   dT    = Tout - DINT(Tout)
%   IF (DABS(Tout-Tref) .lt. eps) THEN
%       QT = QT_i(ivec,1)          <- exact Tref match: return first element
%   ELSE
%       QT = QT_i(ivec,Itemp) + dT*(QT_i(ivec,Itemp+1) - QT_i(ivec,Itemp))
%   IF (QT_i(ivec,Itemp) .lt. 0) QT = -1  <- sentinel for missing data
%
% CRITICAL FIX vs previous version:
%   Previous code used round() for DINT - WRONG
%   FORTRAN DINT truncates toward zero -> MATLAB fix()
%   Example: T=500.7K -> FORTRAN: DINT(500.7)=500, Itemp=205
%                      -> old MATLAB: round(500.7)=501, Itemp=206  WRONG
%                      -> new MATLAB: fix(500.7)=500,   Itemp=205  CORRECT
%
% Inputs:
%   Tgrid - temperature grid vector (used for bounds check only)
%   Qgrid - partition sum row vector for one species (1 x Ntemp)
%   Tval  - temperature at which to evaluate Q (scalar, Kelvin)
%
% Output:
%   QT    - interpolated partition function value (scalar)

    % FORTRAN DATA constants - exact values
    Tref   = 296.0;
    epsTol = 1.0e-5;

    % Ensure Qgrid is a row vector
    Qgrid = Qgrid(:)';
    Ntemp = numel(Qgrid);

    % FORTRAN: IF (DABS(Tout-Tref) .lt. eps) THEN QT = QT_i(ivec,1)
    if abs(Tval - Tref) < epsTol
        QT = Qgrid(1);
        return;
    end

    % FORTRAN: Itemp = DINT(Tout) - 295
    % DINT truncates toward zero -> MATLAB fix()
    % NOT round() - this was the bug in the previous version
    dintT  = fix(Tval);        % DINT(Tout): truncate toward zero
    Itemp  = dintT - 295;      % 1-based index into QT_i columns

    % FORTRAN: dT = Tout - DINT(Tout)
    dT = Tval - double(dintT);

    % Bounds check: if index out of range use linear extrapolation fallback
    % This handles edge temperatures outside the 705-point grid
    if Itemp < 1 || (Itemp + 1) > Ntemp
        QT = interp1(Tgrid(:), Qgrid(:), Tval, 'linear', 'extrap');
        return;
    end

    % FORTRAN sentinel: IF (QT_i(ivec,Itemp) .lt. 0.d0) QT = -1
    if Qgrid(Itemp) < 0.0
        QT = -1.0;
        return;
    end

    % FORTRAN interpolation:
    % QT = QT_i(ivec,Itemp) + dT*(QT_i(ivec,Itemp+1) - QT_i(ivec,Itemp))
    QT = Qgrid(Itemp) + dT * (Qgrid(Itemp + 1) - Qgrid(Itemp));

end