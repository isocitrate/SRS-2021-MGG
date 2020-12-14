function varargout = process_ra_beta_bandpass(varargin)
% R. Avci
% May 2018
% @=======================================================================
%macro_methodcall;
eval(macro_method)
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()%#ok<DEFNU>
        % Description the process
    sProcess.Comment     = 'Filter-Bandpass';
    sProcess.FileTag     = '| Bandpass';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'SARA';
    sProcess.Index       = 1101;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % Definition of the options
    % === Low bound
    sProcess.options.highpass.Comment = 'Lower cutoff frequency:';
    sProcess.options.highpass.Type    = 'value';
    sProcess.options.highpass.Value   = {0.5,'Hz ',[]};
    % === High bound
    sProcess.options.lowpass.Comment = 'Upper cutoff frequency:';
    sProcess.options.lowpass.Type    = 'value';
    sProcess.options.lowpass.Value   = {50,'Hz ',[]};
end

%% ===== GET OPTIONS =====
function [hp, lp] = GetOptions(sProcess)
    hp = sProcess.options.highpass.Value{1};
    lp = sProcess.options.lowpass.Value{1};
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)%#ok<DEFNU>
    %[hp, lp] = GetOptions(sProcess);
    %disp(hp)
    %Comment = ['Bandpass Filter: ' hp '-' lp ' Hz'];
    Comment = 'Bandpass Filter';
end

%% ===== RUN =====
function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    clc;
    clear sOutput;
    sOutput = {};
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
    megID = strcmp(ChannelType, 'MEG');
    megData = Data(megID,:);
    
    % Remove zeros at the end of meg data
    rng = find(sum(megData(:,:)==0,1)<size(megData,1));
    Data = Data(:,rng);               
    % Get Filter Parameters
    [hp, lp] = GetOptions(sProcess);
    % Bandpass Filter
    DataOut = sara_band_pass_filter(hp,lp,Data,sf);
    
    % Create RAW File
    ProtocolMat = bst_get('ProtocolInfo');
    ProtocolDataDir = ProtocolMat.STUDIES;
    [FileDir, FileName] = fileparts(BSFileName);
    [fdir, fName] = fileparts(FileDir);
    SubjectName = fdir;
    StudyName = fName(5:end);
    RawFileName = strcat(ProtocolDataDir, '/', FileDir, '/', ...
                           StudyName, '_Bandpass.raw');
    clear datamat;
    dmat.data=DataOut;
    dmat.sfreq=sf;
    dmat.Comment='data_bandpass';
    [datamat]=CreatedbBrainstorm(dmat);
    % Load Channel Information
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    ChannelMat  = in_bst_channel(ChannelFile);
    
    channelmat = db_template('channelmat');
    channelmat.Comment = ChannelMat.Comment;
    channelmat.Channel = ChannelMat.Channel;
    FileFormat = 'EEG-EGI-RAW';
    out_data_egi_raw(datamat, RawFileName);
    % export to BS 
    OutputFiles = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                      'subjectname',    SubjectName, ...
                      'datafile',       {RawFileName, FileFormat}, ...
                      'channelreplace', 0, ...
                      'channelalign',   0, ...
                      'evtmode',        'value');
    save(strcat(ProtocolDataDir, '/', OutputFiles.ChannelFile), '-struct', 'channelmat')
    bst_process('CallProcess', 'process_import_channel', OutputFiles, [], ...
                'channelfile', {OutputFiles.ChannelFile, 'BST'});
    bst_progress('stop');
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
%% Bandpass Filter
function D = sara_band_pass_filter(cf1,cf2,D,sf)
    if size(D,1) < size(D,2)
        D = D';
    end
    [b,a] = butter(4,[cf1/(sf/2) cf2/(sf/2)]);
    for i = 1 : size(D,2)
        D(:,i) = D(:,i) - mean(D(:,i));
    end
    D = filtfilt(b,a,D);
    if size(D,2) < size(D,1)
        D = D';
    end
end
%% Peak Detection
function dout = rha(din)
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
    [index, ~]=extract_peaks(din,thr1,win);
    [mis1, extra1] = evalRR(index, Fs);
    error = mis1+extra1;
    min_error = error;
    thr = thr1;
    min_thr = thr;
    for k = 1:100
        thr = thr*1.02;    % raise threshold by 2%
        [index, ~]=extract_peaks(din,thr,win);
        [missed, extra] = evalRR(index, Fs);
        error = extra + missed;
        if (error < min_error)
            min_thr = thr;
            min_error = error;
        end
    end
    [index, amp]=extract_peaks(din,min_thr,win);
    % when here, best score is achieved
    damp = (amp-median(amp))/median(amp);
    mark = find(abs(damp) < .5); % keep if amplitude is within 50% of median
    index = index(mark);
    amp = amp(mark);
    %[missed, extra] = evalRR(index, Fs);
    %disp([missed extra])
end
function [yi y]=extract_peaks(din, thr, win)
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
function [missed, extra]=evalRR(yi, Fs)
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
%% Create Raw Data for Brainstorm 
function [datamat]=CreatedbBrainstorm(datamat)
    data=datamat.data;
    sfreq=datamat.sfreq;
    filename=datamat.Comment;
    if isfield(datamat,'Time')
      Time=datamat.Time;
    else
      Time=(0:size(data,2)-1)/sfreq;
    end
    nSensors = size(data,1);
    datamat = db_template('datamat');
    datamat.Comment = filename;
    datamat.Device='CTF';
    datamat.Time = Time;
    datamat.ChannelFlag = ones(nSensors,1);
    datamat.F = zeros(nSensors, length(datamat.Time));
    datamat.F = data;
    datamat.sf = sfreq;
end