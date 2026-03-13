function optDepth = INPUT(optDepth, nMolecules, hitranFiles, NLINES, MOL, waveMin, waveMax, fineGridStep, ilsHalfWidth, ...
                    temperature, pathLength, pressure, nGridPoints, QT_i, Tgrid, ISOVEC, ISONM)
% INPUT - Read reduced HITRAN files, apply T-corrections, build optDepth array
% Exact translation of FORTRAN SUBROUTINE INPUT

    C2   = 1.438786;
    Tref = 296.0;
    TTref  = temperature * Tref;
    DTTref = temperature - Tref;
    WVLW   = waveMin - 45.0;
    WVHG   = waveMax + 45.0;
    IFLG   = 1;

    for K = 1:nMolecules

        fid = fopen(hitranFiles{K}, 'r');
        if fid == -1
            error('INPUT: cannot open HITRAN file: %s', hitranFiles{K});
        end

        lines_read      = 0;
        lines_in_window = 0;
        lines_skipped   = 0;
        lines_added     = 0;

        for I = 1:NLINES(K)

            tline = fgetl(fid);
            if ~ischar(tline), break; end
            lines_read = lines_read + 1;

            if length(tline) < 71
                tline = [tline, repmat(' ', 1, 71 - length(tline))];
            end

            tline = strrep(tline, 'D', 'E');
            tline = strrep(tline, 'd', 'e');

            Molec = str2double(strtrim(tline(1:2)));
            ISO   = str2double(strtrim(tline(3:3)));
            V1    = str2double(strtrim(tline(4:15)));
            S1    = str2double(strtrim(tline(16:25)));
            GM1   = str2double(strtrim(tline(36:40)));
            E1    = str2double(strtrim(tline(46:55)));
            GN1   = str2double(strtrim(tline(56:59)));

            if any(isnan([Molec, ISO, V1, S1, GM1, E1, GN1]))
                lines_skipped = lines_skipped + 1;
                continue;
            end

            if V1 < WVLW || V1 > WVHG
                continue;
            end
            lines_in_window = lines_in_window + 1;

            if Molec < 1 || Molec > numel(ISOVEC)
                lines_skipped = lines_skipped + 1;
                continue;
            end

            niso_mol = ISONM(Molec);
            if ISO < 1 || ISO > niso_mol
                lines_skipped = lines_skipped + 1;
                continue;
            end

            ivec = ISOVEC(Molec) + ISO;

            if ivec < 1 || ivec > size(QT_i, 1)
                lines_skipped = lines_skipped + 1;
                continue;
            end

            QRF_val = QofT_matlab(Tgrid, QT_i(ivec, :), Tref);
            QT_val  = QofT_matlab(Tgrid, QT_i(ivec, :), temperature);

            if QRF_val <= 0.0
                lines_skipped = lines_skipped + 1;
                continue;
            end

            if E1 < 0.0
                lines_skipped = lines_skipped + 1;
                continue;
            end

            QRF_T_val = QRF_val / QT_val;

            Rad = (1.0 - exp(-C2 * V1 / temperature)) / ...
                  (1.0 - exp(-C2 * V1 / Tref));
            lineStrength  = S1 * QRF_T_val * Rad * exp(C2 * E1 * DTTref / TTref);
            halfWidth = GM1 * pressure * (Tref / temperature)^GN1;

            optDepth = OPTD(V1, lineStrength, halfWidth, K, optDepth, nGridPoints, nMolecules, ...
                      waveMin, waveMax, fineGridStep, temperature, pathLength, IFLG);
            IFLG = 2;
            lines_added = lines_added + 1;

        end

        fclose(fid);

        % Diagnostic output per molecule
        fprintf(1, ['INPUT molecule %d (%s): ' ...
            'read=%d  in_window=%d  skipped=%d  added_to_OD=%d\n'], ...
            K, 'mol', lines_read, lines_in_window, ...
            lines_skipped, lines_added);

    end

end