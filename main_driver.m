function main_driver()
% main_driver - FTIR Data Reduction Program
% Exact translation of FORTRAN main driver from QFTIR 2006.FOR
%
% TWO SEPARATE GRID PARAMETERS (critical distinction):
%   fineGridStep = 0.015 cm^-1  -> fine optDepth grid step (INPUT/OPTD/INFR/CONV interior)
%   ilsHalfWidth = 1.28565 cm^-1 -> ILS half-width (APND/CONV only)
%
% nGridPoints = fix((waveMax-waveMin+60)/fineGridStep) = fix(416/0.015) = 27733
% nConvPoints = fix(30/fineGridStep) = 2000 integration points per data point in CONV

    % =========================================================
    % FORTRAN BLOCK DATA ISOVECT
    % =========================================================
    ISOVEC = [ ...
        0,   6,  14,  19,  24,  30,  33, ...
       36,  39,  41,  42,  44,  45,  48, ...
       49,  51,  53,  54,  56,  61,  64, ...
       66,  67,  70,  72,  73,  75,  76, ...
       77,  78,  79,  82,  83,  84,  85, ...
       87,  88,  90,  92];

    ISONM = [ ...
        6,   8,   5,   5,   6,   3,   3, ...
        3,   2,   1,   2,   1,   3,   1, ...
        2,   2,   1,   2,   5,   3,   2, ...
        1,   3,   2,   1,   2,   1,   1, ...
        1,   1,   3,   1,   1,   1,   2, ...
        1,   2,   2,   1];

    % =========================================================
    % GUI SETUP
    % =========================================================
    hFig = figure('Name','FTIR Analysis Setup','NumberTitle','off', ...
        'MenuBar','none','Toolbar','none', ...
        'Position',[300 300 950 700],'Resize','off');
    movegui(hFig,'center');

    uicontrol(hFig,'Style','text','Position',[30 640 160 20], ...
        'String','Path Length (cm):','HorizontalAlignment','right');
    hPathLength = uicontrol(hFig,'Style','edit', ...
        'Position',[200 640 80 25],'String','5');

    uicontrol(hFig,'Style','text','Position',[300 640 160 20], ...
        'String','Pressure (atm):','HorizontalAlignment','right');
    hPressure = uicontrol(hFig,'Style','edit', ...
        'Position',[470 640 80 25],'String','1.0');

    uicontrol(hFig,'Style','text','Position',[560 640 160 20], ...
        'String','Temperature (K):','HorizontalAlignment','right');
    hTemperature = uicontrol(hFig,'Style','edit', ...
        'Position',[730 640 80 25],'String','296.2');

    uicontrol(hFig,'Style','text','Position',[30 600 190 20], ...
        'String','Min Wavenumber (cm^{-1}):','HorizontalAlignment','right');
    hMinWn = uicontrol(hFig,'Style','edit', ...
        'Position',[230 600 80 25],'String','1850');

    uicontrol(hFig,'Style','text','Position',[330 600 190 20], ...
        'String','Max Wavenumber (cm^{-1}):','HorizontalAlignment','right');
    hMaxWn = uicontrol(hFig,'Style','edit', ...
        'Position',[530 600 80 25],'String','2206');

    % ilsHalfWidth = ILS half-width — matches data point spacing
    uicontrol(hFig,'Style','text','Position',[30 560 210 20], ...
        'String','ILS Half-Width (cm^{-1}):','HorizontalAlignment','right');
    hDWV = uicontrol(hFig,'Style','edit', ...
        'Position',[250 560 80 25],'String','1.28565');

    % fineGridStep = fine optDepth grid step — much smaller than data spacing
    uicontrol(hFig,'Style','text','Position',[350 560 210 20], ...
        'String','Fine Grid Step (cm^{-1}):','HorizontalAlignment','right');
    hDVI = uicontrol(hFig,'Style','edit', ...
        'Position',[570 560 80 25],'String','0.015');

    uicontrol(hFig,'Style','text','Position',[30 520 210 20], ...
        'String','Initial Baseline a1:','HorizontalAlignment','right');
    hScale = uicontrol(hFig,'Style','edit', ...
        'Position',[250 520 80 25],'String','1.0');

    uicontrol(hFig,'Style','text','Position',[350 520 210 20], ...
        'String','Initial Baseline b1:','HorizontalAlignment','right');
    hOffset = uicontrol(hFig,'Style','edit', ...
        'Position',[570 520 80 25],'String','-1e-7');

    uicontrol(hFig,'Style','text','Position',[30 480 170 20], ...
        'String','FTIR Data Folder:','HorizontalAlignment','right');
    hFTIRFolder = uicontrol(hFig,'Style','edit', ...
        'Position',[210 480 520 25],'Enable','inactive', ...
        'HorizontalAlignment','left');
    uicontrol(hFig,'Style','pushbutton','Position',[740 480 90 25], ...
        'String','Browse...','Callback',@browseFTIRFolder);

    uicontrol(hFig,'Style','text','Position',[30 440 170 20], ...
        'String','HITRAN Folder:','HorizontalAlignment','right');
    hHITRANFolder = uicontrol(hFig,'Style','edit', ...
        'Position',[210 440 520 25],'Enable','inactive', ...
        'HorizontalAlignment','left');
    uicontrol(hFig,'Style','pushbutton','Position',[740 440 90 25], ...
        'String','Browse...','Callback',@browseHITRANFolder);

    uicontrol(hFig,'Style','text','Position',[30 400 170 20], ...
        'String','Partition Sum List File:','HorizontalAlignment','right');
    hPFListFile = uicontrol(hFig,'Style','edit', ...
        'Position',[210 400 520 25],'Enable','inactive', ...
        'HorizontalAlignment','left');
    uicontrol(hFig,'Style','pushbutton','Position',[740 400 90 25], ...
        'String','Browse...','Callback',@browsePFListFile);

    uicontrol(hFig,'Style','text','Position',[30 360 170 20], ...
        'String','Output Folder:','HorizontalAlignment','right');
    hOutFolder = uicontrol(hFig,'Style','edit', ...
        'Position',[210 360 520 25],'Enable','inactive', ...
        'HorizontalAlignment','left');
    uicontrol(hFig,'Style','pushbutton','Position',[740 360 90 25], ...
        'String','Browse...','Callback',@browseOutFolder);

    uicontrol(hFig,'Style','text','Position',[30 320 170 20], ...
        'String','Summary Output File:','HorizontalAlignment','right');
    hNAFI = uicontrol(hFig,'Style','edit', ...
        'Position',[210 320 300 25],'String','results.O');

    uicontrol(hFig,'Style','text','Position',[30 290 870 20], ...
        'String', ...
        ['Molecules (max 6). Mol ID: 1=H2O 2=CO2 3=O3 4=N2O 5=CO ' ...
         '6=CH4 7=O2 8=NO 9=SO2 10=NO2 | Enter Name | Edit Init PP | ' ...
         'Num Lines auto-filled.'], ...
        'HorizontalAlignment','left');
    hMolTable = uitable(hFig, ...
        'Data',{}, ...
        'ColumnName',{'HITRAN File','Mol ID','Name','Init PP (atm)','Num Lines'}, ...
        'ColumnEditable',[false true true true false], ...
        'ColumnWidth',{220,55,110,110,90}, ...
        'Position',[30 140 710 145]);

    uicontrol(hFig,'Style','pushbutton','Position',[760 250 130 35], ...
        'String','Save Settings','FontWeight','bold', ...
        'Callback',@saveSettingsCallback);
    uicontrol(hFig,'Style','pushbutton','Position',[760 200 130 35], ...
        'String','Proceed','FontWeight','bold', ...
        'Callback',@proceedCallback);
    uicontrol(hFig,'Style','pushbutton','Position',[760 150 130 35], ...
        'String','Cancel','FontWeight','bold', ...
        'Callback',@(~,~)close(hFig));

    % =========================================================
    % GUI CALLBACKS
    % =========================================================
    function browseFTIRFolder(~,~)
        fp = uigetdir('','Select FTIR Data Folder');
        if isequal(fp,0), return; end
        set(hFTIRFolder,'String',fp);
    end

    function browseHITRANFolder(~,~)
        fp = uigetdir('','Select HITRAN Folder');
        if isequal(fp,0), return; end
        set(hHITRANFolder,'String',fp);

        hitran_files = [dir(fullfile(fp,'*.RED')); ...
                        dir(fullfile(fp,'*.dat'))];
        if isempty(hitran_files)
            warndlg('No .RED or .dat files found.','No HITRAN Files');
            return;
        end
        if numel(hitran_files) > 6
            warndlg(sprintf('%d files found. Only first 6 used.', ...
                numel(hitran_files)),'nMolecules Limit');
            hitran_files = hitran_files(1:6);
        end

        NM_gui    = numel(hitran_files);
        tableData = cell(NM_gui, 5);
        for ii = 1:NM_gui
            fname    = hitran_files(ii).name;
            fullpath = fullfile(fp, fname);
            fid_c    = fopen(fullpath,'r');
            nl = 0;
            if fid_c ~= -1
                while true
                    tl = fgetl(fid_c);
                    if ~ischar(tl), break; end
                    if ~isempty(strtrim(tl)), nl = nl+1; end
                end
                fclose(fid_c);
            end
            tableData{ii,1} = fname;
            tableData{ii,2} = 0;
            tableData{ii,3} = '';
            tableData{ii,4} = 0.1;
            tableData{ii,5} = nl;
        end
        set(hMolTable,'Data',tableData);
    end

    function browsePFListFile(~,~)
        [fn,fp] = uigetfile({'*.txt;*.dat','Partition Sum List'}, ...
            'Select Input_file_data.txt');
        if isequal(fn,0), return; end
        set(hPFListFile,'String',fullfile(fp,fn));
    end

    function browseOutFolder(~,~)
        fp = uigetdir('','Select Output Folder');
        if isequal(fp,0), return; end
        set(hOutFolder,'String',fp);
    end

    function saveSettingsCallback(~,~)
        settings.PathLength   = get(hPathLength,  'String');
        settings.Pressure     = get(hPressure,    'String');
        settings.Temperature  = get(hTemperature, 'String');
        settings.MinWn        = get(hMinWn,        'String');
        settings.MaxWn        = get(hMaxWn,        'String');
        settings.ilsHalfWidth          = get(hDWV,          'String');
        settings.fineGridStep          = get(hDVI,          'String');
        settings.BaselineA1   = get(hScale,        'String');
        settings.BaselineB1   = get(hOffset,       'String');
        settings.FTIRFolder   = get(hFTIRFolder,   'String');
        settings.HITRANFolder = get(hHITRANFolder, 'String');
        settings.PFListFile   = get(hPFListFile,   'String');
        settings.OutFolder    = get(hOutFolder,    'String');
        settings.summaryFileName         = get(hNAFI,         'String');
        settings.MolTable     = get(hMolTable,     'Data');
        [fn,fp] = uiputfile({'*.mat','Settings (*.mat)'},'Save Settings');
        if isequal(fn,0), return; end
        save(fullfile(fp,fn),'settings');
        msgbox(['Saved to: ',fullfile(fp,fn)],'Saved');
    end

    function proceedCallback(~,~)
        uiresume(gcbf);
    end

    uiwait(hFig);
    if ~ishandle(hFig)
        error('main_driver: user cancelled.');
    end

    % =========================================================
    % COLLECT AND VALIDATE GUI INPUTS
    % =========================================================
    mol_table_data = get(hMolTable,'Data');
    if isempty(mol_table_data)
        errordlg('No molecules defined.','Input Error'); return;
    end

    for k = 1:size(mol_table_data,1)
        molid = mol_table_data{k,2};
        if isempty(molid)||~isnumeric(molid)||isnan(molid)||molid<1||molid>39
            errordlg(sprintf('Mol ID row %d invalid.',k),'Input Error'); return;
        end
        if isempty(strtrim(mol_table_data{k,3}))
            errordlg(sprintf('Name row %d empty.',k),'Input Error'); return;
        end
        pp = mol_table_data{k,4};
        if isempty(pp)||~isnumeric(pp)||isnan(pp)||pp<=0
            errordlg(sprintf('Init PP row %d invalid.',k),'Input Error'); return;
        end
        nl = mol_table_data{k,5};
        if isempty(nl)||~isnumeric(nl)||isnan(nl)||nl<1
            errordlg(sprintf('Num Lines row %d invalid.',k),'Input Error'); return;
        end
    end

    pathLength   = str2double(get(hPathLength,  'String'));
    pressure  = str2double(get(hPressure,    'String'));
    temperature   = str2double(get(hTemperature, 'String'));
    waveMin = str2double(get(hMinWn,       'String'));
    waveMax = str2double(get(hMaxWn,       'String'));
    ilsHalfWidth  = str2double(get(hDWV,         'String'));  % ILS half-width
    fineGridStep  = str2double(get(hDVI,         'String'));  % fine optDepth grid step
    A1_0 = str2double(get(hScale,       'String'));
    B1_0 = str2double(get(hOffset,      'String'));

    if isnan(pathLength)  || pathLength  <= 0, errordlg('Path Length invalid.','Input Error'); return; end
    if isnan(pressure) || pressure <= 0, errordlg('Pressure invalid.','Input Error'); return; end
    if isnan(temperature)  || temperature  <= 0, errordlg('Temperature invalid.','Input Error'); return; end
    if isnan(waveMin)||isnan(waveMax)||waveMin>=waveMax, errordlg('Wavenumber range invalid.','Input Error'); return; end
    if isnan(ilsHalfWidth) || ilsHalfWidth <= 0, errordlg('ilsHalfWidth invalid.','Input Error'); return; end
    if isnan(fineGridStep) || fineGridStep <= 0, errordlg('fineGridStep invalid.','Input Error'); return; end
    if isnan(A1_0), errordlg('Baseline a1 invalid.','Input Error'); return; end
    if isnan(B1_0), errordlg('Baseline b1 invalid.','Input Error'); return; end

    ftir_folder   = get(hFTIRFolder,   'String');
    hitran_folder = get(hHITRANFolder, 'String');
    pf_list_file  = get(hPFListFile,   'String');
    output_folder = get(hOutFolder,    'String');
    summaryFileName          = get(hNAFI,         'String');

    if ~isfolder(ftir_folder),   errordlg('FTIR folder invalid.','Input Error');   return; end
    if ~isfolder(hitran_folder), errordlg('HITRAN folder invalid.','Input Error'); return; end
    if ~isfile(pf_list_file),    errordlg('PF list file invalid.','Input Error');  return; end
    if ~isfolder(output_folder), errordlg('Output folder invalid.','Input Error'); return; end
    if isempty(strtrim(summaryFileName)),   errordlg('Output filename empty.','Input Error'); return; end

    delete(hFig);

    % =========================================================
    % PARSE MOLECULE TABLE
    % =========================================================
    nMolecules = size(mol_table_data,1);
    if nMolecules > 6, error('main_driver: nMolecules=%d exceeds maximum of 6.',nMolecules); end

    hitranFiles      = cell(nMolecules,1);
    MOL        = zeros(nMolecules,1);
    nHitranLines = zeros(nMolecules,1);
    paramNames        = cell(nMolecules+2,1);
    AY_init    = zeros(nMolecules,1);

    for k = 1:nMolecules
        hitranFiles{k}      = fullfile(hitran_folder, mol_table_data{k,1});
        MOL(k)        = mol_table_data{k,2};
        paramNames{k}        = mol_table_data{k,3};
        AY_init(k)    = mol_table_data{k,4};
        nHitranLines(k) = mol_table_data{k,5};
        if ~isfile(hitranFiles{k})
            error('main_driver: HITRAN file not found: %s',hitranFiles{k});
        end
    end
    paramNames{nMolecules+1} = ' A ';
    paramNames{nMolecules+2} = ' B ';
    nParams = nMolecules+2;

    OAY       = zeros(nParams,1);
    OAY(1:nMolecules) = AY_init;
    OAY(nMolecules+1) = A1_0;
    OAY(nMolecules+2) = B1_0;

    OXL   = pathLength;
    OPT1  = pressure;
    ODVI  = fineGridStep;
    OWVMN = waveMin;
    OWVMX = waveMax;
    ODWV  = ilsHalfWidth;

    % =========================================================
    % nGridPoints — MUST USE FINE GRID fineGridStep (0.015), NOT ilsHalfWidth
    % FORTRAN: nGridPoints = DINT((waveMax-waveMin+60.D0)/fineGridStep)
    % fineGridStep=0.015 -> nGridPoints = fix(416/0.015) = 27733
    % =========================================================
    nGridPoints = fix((waveMax - waveMin + 60.0) / fineGridStep);
    fprintf(1,'nGridPoints = %d  (fineGridStep=%.4f, ilsHalfWidth=%.5f)\n', nGridPoints, fineGridStep, ilsHalfWidth);

    fprintf(1,'Loading partition functions...\n');
    [QT_i, Tgrid] = QTofi(pf_list_file);
    fprintf(1,'Partition functions loaded.\n');

    ftir_files = [dir(fullfile(ftir_folder,'*.1'  )); ...
                  dir(fullfile(ftir_folder,'*.txt' )); ...
                  dir(fullfile(ftir_folder,'*.dat' ))];
    nScanFiles = numel(ftir_files);
    if nScanFiles==0, error('main_driver: no FTIR files found in: %s',ftir_folder); end
    fprintf(1,'Found %d FTIR data file(s) to process.\n',nScanFiles);

    nafi_path = fullfile(output_folder,summaryFileName);
    fidSummary = fopen(nafi_path,'w');
    if fidSummary==-1, error('main_driver: cannot open: %s',nafi_path); end

    % =========================================================
    % SCAN LOOP
    % =========================================================
    for K = 1:nScanFiles

        [~,base_name,~] = fileparts(ftir_files(K).name);
        scanFilePath  = fullfile(ftir_folder,   ftir_files(K).name);
        coeffOutFile = fullfile(output_folder, [base_name,'_coeff.txt']);
        transmOutFile = fullfile(output_folder, [base_name,'_transmittance.txt']);

        fprintf(1,'\nProcessing scan %d of %d: %s\n',K,nScanFiles,scanFilePath);

        fitParams   = OAY;
        fineGridStep  = ODVI;
        waveMin = OWVMN;
        waveMax = OWVMX;
        ilsHalfWidth  = ODWV;
        chiSq = 0.0;

        [measWavenumber, measTransmittance, ~, nDataPoints] = INDAT(scanFilePath, waveMin, waveMax);
        if nDataPoints==0
            warning('main_driver: no data points - skipping %s',scanFilePath);
            continue;
        end
        fprintf(1,'  nDataPoints = %d data points in spectral range.\n',nDataPoints);

        % Build fine-grid optDepth array
        optDepth = zeros(nGridPoints+1, nMolecules);
        optDepth = INPUT(optDepth, nMolecules, hitranFiles, nHitranLines, MOL, waveMin, waveMax, fineGridStep, ilsHalfWidth, ...
                   temperature, OXL, OPT1, nGridPoints, QT_i, Tgrid, ISOVEC, ISONM);

        % =========================================================
        % DIAGNOSTIC BLOCK — remove once fit confirmed correct
        % =========================================================
        fprintf(1,'\n--- DIAGNOSTIC ---\n');
        fprintf(1,'nGridPoints=%d  fineGridStep=%.4f  ilsHalfWidth=%.5f\n', nGridPoints, fineGridStep, ilsHalfWidth);
        fprintf(1,'optDepth grid size: %d x %d\n', size(optDepth,1), size(optDepth,2));
        for kk=1:nMolecules
            fprintf(1,'Mol %d (%s): max optDepth=%.4e  sum optDepth=%.4e  nz=%d\n', ...
                kk, paramNames{kk}, max(optDepth(:,kk)), sum(optDepth(:,kk)), sum(optDepth(:,kk)>0));
        end

        [~, idx_max] = min(measTransmittance);
        test_idx = unique([1, idx_max, nDataPoints]);
        for temperature = 1:numel(test_idx)
            idx = test_idx(temperature);
            St  = CONV(measWavenumber(idx), OAY, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);
            fprintf(1,'CONV at measWavenumber(%d)=%.4f: T=%.6f\n', idx, measWavenumber(idx), St(1));
            for kk=1:nMolecules
                fprintf(1,'  dT/dp_%s = %.4e\n', paramNames{kk}, St(kk+1));
            end
            fprintf(1,'  dT/da1 = %.4e  dT/db1 = %.4e\n', St(nMolecules+2), St(nMolecules+3));
        end

        % Sensitivity check at max absorption point
        fprintf(1,'Sensitivity check at measWavenumber(%d)=%.4f:\n', idx_max, measWavenumber(idx_max));
        Sb = CONV(measWavenumber(idx_max), OAY, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);
        for kk=1:nMolecules
            AYp     = OAY;
            AYp(kk) = OAY(kk)*1.01;
            Sp      = CONV(measWavenumber(idx_max), AYp, optDepth, nGridPoints, nMolecules, waveMin, fineGridStep, ilsHalfWidth);
            dTn     = Sp(1) - Sb(1);
            dTa     = Sb(kk+1)*(AYp(kk)-OAY(kk));
            fprintf(1,'  Mol %d (%s): num=%.4e  anal=%.4e', kk, paramNames{kk}, dTn, dTa);
            if abs(dTn)<1e-12
                fprintf(1,'  WARNING zero sensitivity\n');
            else
                fprintf(1,'  ratio=%.6f\n', dTa/dTn);
            end
        end
        fprintf(1,'--- END DIAGNOSTIC ---\n\n');
        % =========================================================
        % END DIAGNOSTIC BLOCK
        % =========================================================

        TA    = zeros(1000,1);
        covariance = zeros(10,10);
        hessian = zeros(10,10);
        chiSqPrev = 0.0;
        paramsTrial   = zeros(nParams,1);
        gradient   = zeros(nParams,1);
        paramStep     = zeros(nParams,1);
        lambdaDamp = -0.1;
        stepImproved     = 0;

        [fitParams,covariance,hessian,chiSq,lambdaDamp,stepImproved,chiSqPrev,paramsTrial,gradient,paramStep] = ...
            MRQMIN(TA,nDataPoints,fitParams,nParams,covariance,hessian,chiSq,lambdaDamp,stepImproved, ...
                   measWavenumber,measTransmittance,optDepth,nGridPoints,nMolecules,waveMin,fineGridStep,ilsHalfWidth,chiSqPrev,paramsTrial,gradient,paramStep);

        fprintf(1,'Initial chiSq = %.7e\n', chiSq);
        fprintf(1,'Initial fitParams    ='); fprintf(1,' %.4e',fitParams); fprintf(1,'\n\n');

        converged = false;
        nIterations       = 41;

        fprintf(1,'Iter  stepImproved  CHISQ_before    CHISQ_after     lambdaDamp\n');
        fprintf(1,'----  --  --------------  --------------  ----------\n');

        for I = 1:40

            CHISQ_before = chiSq;

            [fitParams,covariance,hessian,chiSq,lambdaDamp,stepImproved,chiSqPrev,paramsTrial,gradient,paramStep] = ...
                MRQMIN(TA,nDataPoints,fitParams,nParams,covariance,hessian,chiSq,lambdaDamp,stepImproved, ...
                       measWavenumber,measTransmittance,optDepth,nGridPoints,nMolecules,waveMin,fineGridStep,ilsHalfWidth,chiSqPrev,paramsTrial,gradient,paramStep);

            fprintf(1,'%4d  %2d  %14.7e  %14.7e  %10.4e\n', ...
                I, stepImproved, CHISQ_before, chiSq, lambdaDamp);

            if stepImproved==0 && I>1
                DCHISQ = abs(CHISQ_before - chiSq);
                PCHISQ = 0.0;
                if chiSq>0, PCHISQ = DCHISQ/chiSq; end
                fprintf(1,'      DCHISQ=%.4e  PCHISQ=%.4e\n',DCHISQ,PCHISQ);

                if (DCHISQ<=1.0e-7) || (PCHISQ<=0.1e-7)
                    fprintf(1,'Converged at iteration %d\n',I);
                    fprintf(1,'Final fitParams ='); fprintf(1,' %.4e',fitParams); fprintf(1,'\n');

                    lambdaDamp = 0.0;
                    [fitParams,covariance,hessian,chiSq,lambdaDamp,stepImproved,chiSqPrev,paramsTrial,gradient,paramStep] = ...
                        MRQMIN(TA,nDataPoints,fitParams,nParams,covariance,hessian,chiSq,lambdaDamp,stepImproved, ...
                               measWavenumber,measTransmittance,optDepth,nGridPoints,nMolecules,waveMin,fineGridStep,ilsHalfWidth,chiSqPrev,paramsTrial,gradient,paramStep);

                    nIterations       = I;
                    converged = true;
                    fprintf(1,'%s\n',scanFilePath);
                    OUTPT(measWavenumber,measTransmittance,fitParams,covariance,chiSq,coeffOutFile,transmOutFile,'', ...
                          nMolecules,optDepth,nGridPoints,waveMin,fineGridStep,ilsHalfWidth,paramNames,nIterations,measWavenumber,measTransmittance,nDataPoints,nParams);
                    fprintf(1,'%s\n',scanFilePath);
                    break;
                end
            end

        end

        if ~converged
            fprintf(1,'PROGRAM DID NOT CONVERGE IN 41 ITERATIONS\n');
            fprintf(1,'Final fitParams ='); fprintf(1,' %.4e',fitParams); fprintf(1,'\n');
            nIterations = 41;
            OUTPT(measWavenumber,measTransmittance,fitParams,covariance,chiSq,coeffOutFile,transmOutFile,'', ...
                  nMolecules,optDepth,nGridPoints,waveMin,fineGridStep,ilsHalfWidth,paramNames,nIterations,measWavenumber,measTransmittance,nDataPoints,nParams);
        end

        % Summary file output
        paramsPercent = fitParams * 100.0;
        if converged
            switch nMolecules
                case 1, fprintf(fidSummary,' %2d     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)));
                case 2, fprintf(fidSummary,' %2d     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)));
                case 3, fprintf(fidSummary,' %2d     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)));
                case 4, fprintf(fidSummary,' %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)));
                case 5, fprintf(fidSummary,' %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)), ...
                          paramNames{5},fortran_E10_4(paramsPercent(5)));
                case 6, fprintf(fidSummary,' %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)), ...
                          paramNames{5},fortran_E10_4(paramsPercent(5)), ...
                          paramNames{6},fortran_E10_4(paramsPercent(6)));
            end
        else
            switch nMolecules
                case 1, fprintf(fidSummary,'     %2d     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)));
                case 2, fprintf(fidSummary,'     %2d     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)));
                case 3, fprintf(fidSummary,'     %2d     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)));
                case 4, fprintf(fidSummary,'     %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)));
                case 5, fprintf(fidSummary,'     %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)), ...
                          paramNames{5},fortran_E10_4(paramsPercent(5)));
                case 6, fprintf(fidSummary,'     %2d     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s     %-5s %s\n', ...
                        K,paramNames{1},fortran_E10_4(paramsPercent(1)), ...
                          paramNames{2},fortran_E10_4(paramsPercent(2)), ...
                          paramNames{3},fortran_E10_4(paramsPercent(3)), ...
                          paramNames{4},fortran_E10_4(paramsPercent(4)), ...
                          paramNames{5},fortran_E10_4(paramsPercent(5)), ...
                          paramNames{6},fortran_E10_4(paramsPercent(6)));
            end
        end

    end % scan loop

    fclose(fidSummary);
    fprintf(1,'\n--- All scans processed. Summary: %s ---\n',nafi_path);

end % main_driver


function s = fortran_E10_4(val)
    if val==0.0, s=' 0.0000E+00'; return; end
    sgn=' ';
    if val<0, sgn='-'; val=abs(val); end
    exp10=floor(log10(val));
    mant=val/10^exp10/10.0;
    exp10=exp10+1;
    mant=round(mant*1e4)/1e4;
    if mant>=1.0, mant=mant/10.0; exp10=exp10+1; end
    if exp10>=0
        s=sprintf('%s0.%04dE+%02d',sgn,round(mant*1e4),exp10);
    else
        s=sprintf('%s0.%04dE-%02d',sgn,round(mant*1e4),abs(exp10));
    end
    while length(s)<10, s=[' ',s]; end %#ok<AGROW>
end