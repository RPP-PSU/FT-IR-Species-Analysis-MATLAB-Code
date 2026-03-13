function [QT_i, Tgrid] = QTofi(pf_list_file)
% QTofi - Load partition sum data for all 93 species
% Exact translation of FORTRAN SUBROUTINE QTofi
%
% FORTRAN source:
%   OPEN(2, FILE='C:\QUANT\Partition-Sums\Input_File_data.txt')
%   DO JJ=1,93: READ(2,*) FILDATA(JJ)   <- reads 93 filenames
%   DO ID=1,93:
%     OPEN(30, FILE=FILDATA(ID))
%     DO i=1,705: read(30,*) QT_i(ID,i) <- 705 values per file
%
% DISK LAYOUT (confirmed from user folder screenshot):
%
%   C:\FTIR_DATA\QUANT\                              <- parent_dir
%   C:\FTIR_DATA\QUANT\Partition-Sums\               <- list_folder
%   C:\FTIR_DATA\QUANT\Partition-Sums\Input_file_data.txt
%   C:\FTIR_DATA\QUANT\Partition-Sums\H2O\H2O_161_red.dat
%   C:\FTIR_DATA\QUANT\Partition-Sums\N2O\N2O_446_red.dat
%   ...etc
%
% ENTRIES IN Input_file_data.txt look like:
%   'Partition-Sums/H2O/H2O_161_red.dat'
%
% PATH RESOLUTION:
%   list_folder = fileparts(pf_list_file)
%               = C:\FTIR_DATA\QUANT\Partition-Sums
%
%   parent_dir  = fileparts(list_folder)
%               = C:\FTIR_DATA\QUANT
%
%   fullfile(parent_dir, 'Partition-Sums/H2O/H2O_161_red.dat')
%               = C:\FTIR_DATA\QUANT\Partition-Sums\H2O\H2O_161_red.dat  CORRECT
%
% This matches the FORTRAN original which resolved paths relative to
% C:\QUANT\ (one level above the Partition-Sums folder).
%
% THREE-STRATEGY FALLBACK:
%   1. parent_dir  + relative path  (standard layout - matches your disk)
%   2. list_folder + relative path  (non-standard: list file at root level)
%   3. relative path as-is          (absolute path or already resolvable)
%
% KEY IMPLEMENTATION NOTES:
%   - fscanf() used for reading to match FORTRAN free-format READ(*)
%     which handles variable leading whitespace in .dat files
%   - Exactly 93 files x 705 values enforced (FORTRAN QT_i(93,705))
%   - Single/double quotes stripped from filenames
%   - Forward and back slashes normalized to OS separator
%   - Temperature grid: 296:1000 K (705 points)
%     matching FORTRAN index: Itemp = DINT(Tout)-295, index 1 = 296K
%
% Inputs:
%   pf_list_file - full path to Input_file_data.txt
%                  (selected via GUI browse button)
%
% Outputs:
%   QT_i  - partition sum matrix (93 x 705)
%            rows = species index (1-93 matching Input_file_data.txt order)
%            cols = temperature index (1=296K, 705=1000K)
%   Tgrid - temperature grid vector (705 x 1), values 296:1000 K

    % FORTRAN array dimensions - exact values
    NSPECIES = 93;     % FORTRAN: QT_i(93,705) - 93 species
    NTEMP    = 705;    % FORTRAN: 705 temperature points per species

    % Temperature grid: 1K increments starting at 296K
    % FORTRAN index mapping: Itemp = DINT(Tout)-295
    %   Itemp=1  -> T=296K
    %   Itemp=705 -> T=1000K
    Tgrid = (296 : 296 + NTEMP - 1)';    % 705x1 column vector, 296 to 1000 K

    % =========================================================
    % DETERMINE BASE DIRECTORY FOR RELATIVE PATH RESOLUTION
    %
    % Input_file_data.txt is inside Partition-Sums\
    % Entries in the file start with 'Partition-Sums/...'
    % Therefore base must be ONE LEVEL ABOVE the list file's folder
    % =========================================================
    list_folder = fileparts(pf_list_file);   % folder containing list file
    parent_dir  = fileparts(list_folder);    % one level up = correct base

    % =========================================================
    % OPEN AND READ Input_file_data.txt
    % FORTRAN: OPEN(2,FILE=...) then DO JJ=1,93: READ(2,*) FILDATA(JJ)
    % =========================================================
    fid = fopen(pf_list_file, 'r');
    if fid == -1
        error('QTofi: cannot open partition sum list file:\n  %s', ...
            pf_list_file);
    end

    fileList = cell(NSPECIES, 1);

    for JJ = 1:NSPECIES

        tline = fgetl(fid);

        % Check for premature end of file
        if ~ischar(tline)
            fclose(fid);
            error(['QTofi: unexpected end of file at line %d.\n' ...
                   'Expected %d filenames in: %s'], ...
                JJ, NSPECIES, pf_list_file);
        end

        % Strip leading and trailing whitespace
        tline = strtrim(tline);

        % FORTRAN free-format READ(*) strips surrounding single quotes
        % Entries look like: 'Partition-Sums/H2O/H2O_161_red.dat'
        if length(tline) >= 2
            if (tline(1) == '''' && tline(end) == '''') || ...
               (tline(1) == '"'  && tline(end) == '"')
                tline = tline(2:end-1);
            end
        end

        % Normalize all slashes to OS-native separator
        % Handles both Unix-style (/) and Windows-style (\) in the file
        tline = strrep(tline, '/', filesep);
        tline = strrep(tline, '\', filesep);

        % =========================================================
        % THREE-STRATEGY PATH RESOLUTION
        %
        % Strategy 1 (primary): parent_dir + relative path
        %   C:\FTIR_DATA\QUANT + Partition-Sums\H2O\H2O_161_red.dat
        %   = C:\FTIR_DATA\QUANT\Partition-Sums\H2O\H2O_161_red.dat
        %   This is correct for the confirmed disk layout.
        %
        % Strategy 2 (fallback A): list_folder + relative path
        %   C:\FTIR_DATA\QUANT\Partition-Sums + Partition-Sums\H2O\...
        %   = C:\FTIR_DATA\QUANT\Partition-Sums\Partition-Sums\H2O\...
        %   This would be wrong for standard layout but handles edge cases
        %   where the list file sits at the root of the relative paths.
        %
        % Strategy 3 (fallback B): treat tline as absolute or CWD-relative
        %   Handles cases where Input_file_data.txt contains absolute paths
        % =========================================================
        full_path_parent = fullfile(parent_dir,  tline);
        full_path_same   = fullfile(list_folder, tline);

        if isfile(full_path_parent)
            % Standard layout confirmed - use parent-based path
            fileList{JJ} = full_path_parent;

        elseif isfile(full_path_same)
            % Non-standard layout - list file at root of relative paths
            fileList{JJ} = full_path_same;

        elseif isfile(tline)
            % Absolute path or current working directory relative
            fileList{JJ} = tline;

        else
            % All strategies failed - provide detailed diagnostic message
            fclose(fid);
            error(['QTofi: cannot locate partition sum file for entry %d.\n\n' ...
                   'List file entry : %s\n\n' ...
                   'Tried path 1    : %s\n' ...
                   '  (parent_dir + relative, expected for your layout)\n\n' ...
                   'Tried path 2    : %s\n' ...
                   '  (list_folder + relative)\n\n' ...
                   'Tried path 3    : %s\n' ...
                   '  (as-is)\n\n' ...
                   'Your list file  : %s\n' ...
                   'list_folder     : %s\n' ...
                   'parent_dir      : %s\n\n' ...
                   'Check that the molecule subfolders exist inside:\n' ...
                   '  %s'], ...
                JJ, tline, ...
                full_path_parent, ...
                full_path_same, ...
                tline, ...
                pf_list_file, ...
                list_folder, ...
                parent_dir, ...
                list_folder);
        end

    end

    fclose(fid);

    % =========================================================
    % LOAD PARTITION SUM DATA FROM EACH .dat FILE
    % FORTRAN:
    %   DO ID=1,93
    %     OPEN(30,FILE=FILDATA(ID),STATUS='UNKNOWN')
    %     DO i=1,705
    %       read(30,*) QT_i(ID,i)   <- free-format: handles leading spaces
    %     END DO
    %   END DO
    %
    % FORTRAN free-format READ(30,*) handles variable leading whitespace.
    % Confirmed .dat file format (from sample provided):
    %   "           173.633000"   <- variable leading spaces
    %   "          174.515603"
    % MATLAB load() CANNOT handle this -> must use fscanf()
    % =========================================================
    QT_i = zeros(NSPECIES, NTEMP);

    for ID = 1:NSPECIES

        fid2 = fopen(fileList{ID}, 'r');
        if fid2 == -1
            error('QTofi: cannot open partition sum file:\n  %s', ...
                fileList{ID});
        end

        % fscanf reads all whitespace-delimited floats
        % exactly replicating FORTRAN free-format READ(30,*) behavior
        data = fscanf(fid2, '%f');
        fclose(fid2);

        % Enforce exactly NTEMP=705 values per file
        % matching FORTRAN array declaration QT_i(93,705)
        if numel(data) < NTEMP
            error(['QTofi: file has too few values.\n' ...
                   'File   : %s\n' ...
                   'Found  : %d values\n' ...
                   'Expected: %d values (FORTRAN NTEMP=705)'], ...
                fileList{ID}, numel(data), NTEMP);
        end

        % Store exactly NTEMP values as row ID of QT_i
        % FORTRAN: QT_i(ID,1) through QT_i(ID,705)
        QT_i(ID, :) = data(1:NTEMP)';

    end

    fprintf(1, 'QTofi: loaded %d species x %d temperature points.\n', ...
        NSPECIES, NTEMP);

end
