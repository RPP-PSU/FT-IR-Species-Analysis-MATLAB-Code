function [A, B] = GAUSSJ(A, N, B, M)
% GAUSSJ - Gauss-Jordan elimination with full pivoting
% Exact translation of FORTRAN SUBROUTINE GAUSSJ
%
% FORTRAN source uses:
%   - Full column and row pivoting
%   - pause 'singular matrix' if pivot is zero or if ipiv>1
%   - In-place matrix inversion: A becomes A^-1
%   - B is solved in-place: B becomes A^-1 * B
%
% Inputs:
%   A - N x N matrix (will be overwritten with its inverse)
%   N - dimension of A
%   B - N x M right-hand side matrix
%   M - number of right-hand sides
%
% Outputs:
%   A - N x N inverse of input A
%   B - N x M solution matrix

    NMAX = 50;
    indxc = zeros(1, NMAX);
    indxr = zeros(1, NMAX);
    ipiv  = zeros(1, NMAX);

    % FORTRAN DO 22 i=1,n (main elimination loop)
    for i = 1:N
        big  = 0.0;
        irow = 0;
        icol = 0;

        % Search for pivot: DO 13 j=1,n / DO 12 k=1,n
        for j = 1:N
            if ipiv(j) ~= 1
                for k = 1:N
                    if ipiv(k) == 0
                        if abs(A(j,k)) >= big
                            big  = abs(A(j,k));
                            irow = j;
                            icol = k;
                        end
                    elseif ipiv(k) > 1
                        error('GAUSSJ: singular matrix (ipiv > 1)');
                    end
                end
            end
        end

        ipiv(icol) = ipiv(icol) + 1;

        % Swap rows irow and icol if needed (DO 14, DO 15)
        if irow ~= icol
            for l = 1:N
                dum        = A(irow, l);
                A(irow, l) = A(icol, l);
                A(icol, l) = dum;
            end
            for l = 1:M
                dum        = B(irow, l);
                B(irow, l) = B(icol, l);
                B(icol, l) = dum;
            end
        end

        indxr(i) = irow;
        indxc(i) = icol;

        if A(icol, icol) == 0.0
            error('GAUSSJ: singular matrix (zero pivot)');
        end

        pivinv        = 1.0 / A(icol, icol);
        A(icol, icol) = 1.0;

        % Scale pivot row (DO 16, DO 17)
        for l = 1:N
            A(icol, l) = A(icol, l) * pivinv;
        end
        for l = 1:M
            B(icol, l) = B(icol, l) * pivinv;
        end

        % Eliminate column (DO 21 ll=1,n)
        for ll = 1:N
            if ll ~= icol
                dum        = A(ll, icol);
                A(ll,icol) = 0.0;
                for l = 1:N
                    A(ll, l) = A(ll, l) - A(icol, l) * dum;
                end
                for l = 1:M
                    B(ll, l) = B(ll, l) - B(icol, l) * dum;
                end
            end
        end
    end

    % Unscramble columns in reverse order (DO 24 l=n,1,-1)
    for l = N:-1:1
        if indxr(l) ~= indxc(l)
            for k = 1:N
                dum              = A(k, indxr(l));
                A(k, indxr(l))   = A(k, indxc(l));
                A(k, indxc(l))   = dum;
            end
        end
    end

end
