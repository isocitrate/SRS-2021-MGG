function varargout = process_bl_visualisation(varargin)

eval(macro_method);
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()%#ok<DEFNU>
        % Description the process
    sProcess.Comment     = 'Visualise SQUID signals';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'MGG';
    sProcess.Index       = 1003;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    
    % Definition of the options
    % Options: Time window
    sProcess.options.timewindow.Comment = 'Time interval (0 = First timepoint):';
    sProcess.options.timewindow.Type    = 'range';
    sProcess.options.timewindow.Value   = {[0,30], 'sec', 3} ;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)%#ok<DEFNU>

    Comment = 'Show MGG signals for each SQUID sensor';
end

function [first_timept, second_timept] = GetOptions(sProcess)
    first_timept = sProcess.options.timewindow.Value{1}(1);
    second_timept  = sProcess.options.timewindow.Value{1}(2);
    if second_timept == first_timept || first_timept < 0 || second_timept < 0
        bst_error('Invalid time interval selected', 'Select time interval', 0);
        return;
    end
    if first_timept > second_timept
        temp_first = first_timept;
        first_timept = second_timept;
        second_timept = first_timept;
    end 
end

%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
    clc;
    % Get Data from BS
    BSFileName = sInput.FileName;
    [sStudy, ~, ~, ~] = bst_get('AnyFile', BSFileName);
    % Read BS file for all channels
    DataMat = in_bst_data(BSFileName);
    DataFile = DataMat.F;
    % SamplesBounds = DataFile.prop.samples;
    
    % Note, no SamplesBounds exists here because the sFile created during
    % raw file creation doesn't contain any epochs / events (time points
    % correlated to specific events)
    SamplesBounds = [];
 
    % Load channel file
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    ChannelMat  = in_bst_channel(ChannelFile);
    cFlag = DataMat.ChannelFlag';
    cName = {ChannelMat.Channel.Name};
    cType = {ChannelMat.Channel.Type};

    gradID = find((strcmp(cType, 'Grad')==1 & cFlag==1)==1);
    % Get MEG Data
    [gradData] = getDataFromDB(sProcess, DataFile, SamplesBounds, ...
        ChannelMat, gradID);
    subplot_pos = [50 88 95 84 73 46 19 12 5 16 27 52 69 67 48 31 33 54 81];

    for i = 1:length(cName)
        if cFlag(i) ~= 1
            subplot_pos(i) = [];
        end
    end
    % Get Bad Segments
    BadSegments = panel_record('GetBadSegments', DataFile);                

% Loading ChannelMat data
% SenIdsUsed = [33 3 35 9 15 25 23 26 30 36 5 7 16 17 12 14 18 22 19];
%% Generate synthetic data by finding Euclidean distance and generating time delay

% Initialising
[first_timept, second_timept] = GetOptions(sProcess);

for i = 1:length(gradID)
    grad_sensorName{i} = cName{gradID(i)};
end
timeValue_first = first_timept + DataMat.Time(1);
timeValue_second = second_timept + DataMat.Time(1);
firstIndex_array = abs(DataMat.Time - timeValue_first);
[closest_min, firstIndex] = min(firstIndex_array);
[closest_max, lastIndex] = min(abs(DataMat.Time - timeValue_second));
timeWindow = DataMat.Time(firstIndex : lastIndex);
MGG_sig_in_window = gradData(:, firstIndex:lastIndex);
y_min = min(min(MGG_sig_in_window, [], 2)); 
y_max = max(max(MGG_sig_in_window, [], 2)); 
% [first_timept + DataMat.Time(1) , second_timept + DataMat.Time(1)]
% Scan between those indices for all data values (ie find index)). 
% max(max(gradData, [], 2));    
gradData = gradData*10E12;
for i = 1:length(grad_sensorName)
    sp = subplot(11, 9, subplot_pos(i));
    plot(DataMat.Time, gradData(i, :))
    xlim([first_timept + DataMat.Time(1) , second_timept + DataMat.Time(1)]);
    %ylim([y_min, y_max]);
    str = sprintf('Sensor %s', grad_sensorName{i});
    title(str)
    xlabel('Time (s)');
    xticks ([0 20 40 60]);
    set(gca,'Units','normalized')
    yl = ylabel('MF Magnitude (pT)');    
    % Adjust plot size
    d = get(sp,'Position'); % left bottom width height]
    d(3)=d(3)*1.25; % Scale width up by 3
    d(4)=d(4)*1.75; % Scale height up by 3
    set(sp,'Position',d);
    
     yl_pos = get(yl , 'position');
     yl_pos(1) = -7;
     set(yl , 'position' , yl_pos);
end

    %% ===== Read Data from BS Database =====
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
    %% ===== Save File =====
    function SaveFile(iStudy, BSFileName, DataType, RowNames, TFpsd, ...
            OPTIONS, FreqBands, SurfaceFile, GridLoc, GridAtlas, nAvgFile)
        OutputFiles = {};
        % Create file structure
        FileMat = db_template('timefreqmat');
        FileMat.Comment   = OPTIONS.Comment;
        FileMat.DataType  = DataType;
        FileMat.TF        = TFpsd;
        FileMat.Time      = OPTIONS.TimeVector;
        FileMat.TimeBands = [];
        FileMat.Freqs     = FreqBands;
        FileMat.RowNames  = RowNames;
        FileMat.Measure   = OPTIONS.Measure;
        FileMat.Method    = OPTIONS.Method;
        FileMat.nAvg      = nAvgFile;
        FileMat.SurfaceFile = SurfaceFile;
        FileMat.GridLoc     = GridLoc;
        FileMat.GridAtlas   = GridAtlas;
        FileMat.DataFile = file_short(BSFileName);
        FileMat.Options.Method          = OPTIONS.Method;
        FileMat.Options.Measure         = OPTIONS.Measure;
        FileMat.Options.Output          = OPTIONS.Output;
        
        % Save the file
        if ~isempty(iStudy)
            % Get output study
            sTargetStudy = bst_get('Study', iStudy);
            % Output filename
            fileBase = 'timefreq';
            fileBase = [fileBase '_' OPTIONS.Method];
            FileName = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), fileBase);
            % Save file
            bst_save(FileName, FileMat, 'v6');
            % Add file to database structure
            db_add_data(iStudy, FileName, FileMat);
            % Return new filename
            sOutput{end+1} = FileName;
        % Returns the contents of the file instead of saving them
        else
            sOutput{end+1} = FileMat;
        end
    end
end
