function varargout = process_ra_beta_detect_R_markers(varargin)
% R. Avci
% May 2018
% @=======================================================================
%macro_methodcall;
eval(macro_method)
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()%#ok<DEFNU>
        % Description the process
    sProcess.Comment     = 'Removal-Detect R Markers';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'SARA';
    sProcess.Index       = 1201;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % Definition of the options
    sProcess.options.mMCG.Comment = 'Maternal R Markers';
    sProcess.options.mMCG.Type    = 'checkbox';
    sProcess.options.mMCG.Value   = 0;
    sProcess.options.fMCG.Comment = 'Fetal R Markers';
    sProcess.options.fMCG.Type    = 'checkbox';
    sProcess.options.fMCG.Value   = 0;
end

%% ===== GET OPTIONS =====
function [ismMCG, isfMCG] = GetOptions(sProcess)
    ismMCG = sProcess.options.mMCG.Value;
    isfMCG = sProcess.options.fMCG.Value;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)%#ok<DEFNU>
    Comment = 'R Marker Detection';
end

%% ===== RUN =====
function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    clc;
    clear sOutput;
    sOutput = {};
    % Get process options 
    [ismMCG, isfMCG] = GetOptions(sProcess);
    if ~ismMCG && ~isfMCG
        disp('Terminating the process since mMCG/fMCG Selection was not done!')
        return
    elseif ismMCG && isfMCG
        disp('Terminating the process since both mMCG and fMCG were selected!')
        return
    end
    % Get Data from BS
    BSFileName = sInputs.FileName;
    [sStudy, iStudy, ~, DataType] = bst_get('AnyFile', BSFileName);
    % Read BS file for all channels
    DataMat = in_bst_data(BSFileName);
    DataFile = DataMat.F;
    sf = 1 ./ (DataMat.Time(2) - DataMat.Time(1));
    SamplesBounds = DataFile.prop.samples;
    % Load channel file
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    ChannelMat  = in_bst_channel(ChannelFile);
    ChannelFlag = DataMat.ChannelFlag';
    ChannelName = {ChannelMat.Channel.Name};
    ChannelType = {ChannelMat.Channel.Type};
    Data = getDataFromDB(sProcess, DataFile, SamplesBounds, ...
                              ChannelMat, 1:length(ChannelFlag));
    disp(size(Data))
    megID = strcmp(ChannelType, 'MEG');
    megData = Data(megID,:);
    % Compute R markers
    if ismMCG
        din = sum(sara_rha(megData),1);
        iRm = sara_find_peaks_M(din, sf);
        if size(iRm,1)>size(iRm,2)
            iRm = iRm';
        end
        events(1,1).label = 'mMCG';
        events(1,1).color = [0 0 1];
        events(1,1).epochs = ones(size(iRm));
        events(1,1).samples = iRm;
        events(1,1).times = (iRm-1)/sf;
        events(1,1).reactTimes = [];
        events(1,1).select = 1;
    elseif isfMCG
        din = sum(sara_rha(megData),1);
        iRf = sara_find_peaks_F(din, sf);
        if size(iRf,1)>size(iRf,2)
            iRf = iRf';
        end
        events(1,1).label = 'fMCG';
        events(1,1).color = [0 1 0];
        events(1,1).epochs = ones(size(iRf));
        events(1,1).samples = iRf;
        events(1,1).times = (iRf-1)/sf;
        events(1,1).reactTimes = [];
        events(1,1).select = 1;
    end
    DataFile.events = events;
    DataMat.F = DataFile;
    % Save file definition
    bst_save(file_fullpath(BSFileName), DataMat, 'v6', 1);
end
%% Functions
%% Get Data from BS database
function [data, time] = getDataFromDB(sProcess, sFile, SamplesBounds, ...
                                      ChannelMat, ChID)
    % Reading options
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ImportMode = 'Time';
    ImportOptions.Resample   = 0;
    ImportOptions.UseCtfComp = 1;
    ImportOptions.UseSsp     = 1;
    ImportOptions.RemoveBaseline = 'no';
    ImportOptions.DisplayMessages = 0;

    % Read data
    [data, time] = in_fread(sFile, ChannelMat, 1, ...
                            SamplesBounds, ChID, ...
                            ImportOptions);
end
%% Peak Detection
function dout = sara_rha(din)
    % JD Wilson 8-25-2008
    dout = zeros(size(din)); 
    dout(1:end-1,:) = abs(diff(hilbert(din)));
end
function [index, amp] = sara_find_peaks_M(din,Fs)
    % JD Wilson 8-25-2008
    win = 80;
    din(end-win-2:end) = din(end-win-2);      % fix last few data points
    din(1:win+2) = din(win+2);           % fix first few data points
    [nh, xh]=hist(din);
    [~, ih] = max(nh);
    thr1 = xh(ih);
    [index, ~]=sara_extract_peaks(din,thr1,win);
    [mis1, extra1] = sara_evalRR(index, Fs);
    error = mis1+extra1;
    min_error = error;
    thr = thr1;
    min_thr = thr;
    for k = 1:100
        thr = thr*1.02;    % raise threshold by 2%
        [index, ~]=sara_extract_peaks(din,thr,win);
        [missed, extra] = sara_evalRR(index, Fs);
        error = extra + missed;
        if (error < min_error)
            min_thr = thr;
            min_error = error;
        end
    end
    [index, amp]=sara_extract_peaks(din,min_thr,win);
    % when here, best score is achieved
    damp = (amp-median(amp))/median(amp);
    mark = find(abs(damp) < .5); % keep if amplitude is within 50% of median
    index = index(mark);
    amp = amp(mark);
    %[missed, extra] = evalRR(index, Fs);
    %disp([missed extra])
end
function [index, amp] = sara_find_peaks_F(din, Fs)
    % JD Wilson 8-25-2008
    win = 20;               % window for peak search
    din(end-win-2:end) = din(end-win-2);      % fix last few data points
    din(1:win+2) = din(win+2);           % fix first few data points
    [nh, xh]=hist(din);
    [~, ih] = max(nh);
    thr1 = xh(ih);
    [index ~]=sara_extract_peaks(din,thr1,win);
    [mis1 extra1] = sara_evalRR(index, Fs);
    error = mis1+extra1;
    min_error = error;
    thr = thr1;
    min_thr = thr;
    for n = 1:85   %85
        thr = thr*1.02;    % raise threshold by 2%
        [index ~]=sara_extract_peaks(din,thr,win);
        [missed, extra] = sara_evalRR(index, Fs);
        error = extra + missed;
        %disp([n extra missed])
        if (error < min_error)
            min_thr = thr;
            min_error = error;
            %disp([n extra missed])
        end
        if error == 0
            break;
        end
    end
    [index amp]=sara_extract_peaks(din,min_thr,win);
    % when here, best score is achieved
    % amplitude filter
    damp = (amp-median(amp))/median(amp);
    markHi = damp < 2;
    markLo = damp > -0.5;
    mark = find(markHi & markLo);
    index = index(mark);
    amp = amp(mark);
end
function [yi y]=sara_extract_peaks(din, thr, win)
    mark = din>thr;
    mark = find(not(logical(mark(1:end-1))) & logical(mark(2:end)));
    nm = length(mark);
    size(din);
    yi = zeros(nm,1);   
    y = zeros(nm,1);
    for n=1:nm
        p = mark(n);
        [a b]=max(din(p-win:p+win));
        yi(n) = b+p-win-1;      % use -1 to shift left one time position
        y(n)=a;
    end
    while(1)    % remove duplicate peaks
        k = find(diff(yi) == 0);
        if isempty(k)
            break
        end
        y(k)=[]; yi(k)=[];
    end
end
function [missed, extra]=sara_evalRR(yi, Fs)
    dyi = diff(yi)/Fs;
    mrr = median(dyi);
        % See Hilbert paper for values
    if ((mrr < 0.4) || (mrr > 1.2))
        mrr = .69;    % default if out of range
    end
    k = find(dyi>1.2);
    M = length(k);
    missed = (sum(dyi(k))/mrr) - M;
    k = find(dyi<0.4);
    E = length(k);
    extra = E - (sum(dyi(k)/mrr));
    %disp([mrr  missed  extra])
end
