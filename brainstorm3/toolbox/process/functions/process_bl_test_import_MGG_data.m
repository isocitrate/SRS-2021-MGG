function varargout = process_bl_test_import_MGG_data( varargin )

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Read RAW data';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'MGG';
    sProcess.Index       = 1001;
    sProcess.isSeparator = 1;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    % Option: Subject name
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = 'NewSubject';
    % Option: File to import
    sProcess.options.datafile.Comment = 'File to import:';
    sProcess.options.datafile.Type    = 'datafile';
    sProcess.options.datafile.Value   = {...
        '', ...                                % Filename
        '', ...                                % FileFormat
        'open', ...                            % Dialog type: {open,save}
        'Open raw EEG/MEG recordings...', ...  % Window title
        'ImportData', ...                      % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                          % Selection mode: {single,multiple}
        'files_and_dirs', ...                  % Selection mode: {files,dirs,files_and_dirs}
        bst_get('FileFilters', 'raw'), ...    % Get all the available file formats
        'DataIn'};                             % DefaultFormats
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    % Get subject name
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'Subject name is empty.');
        return
    end
    % Get filename to import
    FileName   = sProcess.options.datafile.Value{1};
    FileFormat = sProcess.options.datafile.Value{2};
    if isempty(FileName)
        bst_report('Error', sProcess, [], 'No file selected.');
        return
    end    
    % ===== IMPORT FILES =====
    % Get subject 
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Create subject is it does not exist yet
    if isempty(sSubject)
        [sSubject, iSubject] = db_add_subject(SubjectName);
    end
    % Import options
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ChannelReplace  = 0;
    ImportOptions.ChannelAlign    = 0;
    ImportOptions.DisplayMessages = 0;
    ImportOptions.EventsMode      = 'ignore';
    ImportOptions.EventsTrackMode = 'value';
    % Import link to raw
    OutputFiles{1} = ra_test_import_raw(FileName, FileFormat, iSubject, ImportOptions);
end

function OutputDataFile = ra_test_import_raw(RawFiles, FileFormat, iSubject, ImportOptions)
% IMPORT_RAW: Create a link to a raw file in the Brainstorm database.
%
% USAGE:  OutputDataFile = import_raw(RawFiles, FileFormat, iSubject, ImportOptions)
%         OutputDataFile = import_raw(RawFiles, FileFormat, iSubject)
%         OutputDataFile = import_raw(RawFiles, FileFormat)
%         OutputDataFile = import_raw()
%
% INPUTS:
%     - RawFiles      : Full path to the file to import in database
%     - FileFormat    : String representing the file format (CTF, FIF, 4D, ...)
%     - iSubject      : Subject indice in which to import the raw file
%     - ImportOptions : Structure that describes how to import the recordings
%       => Fields used: ChannelAlign, ChannelReplace, DisplayMessages, EventsMode, EventsTrackMode

%% ===== PARSE INPUT =====
if (nargin < 4) || isempty(ImportOptions)
    ImportOptions = db_template('ImportOptions');
end
if (nargin < 3)
    iSubject = [];
end
if (nargin < 2)
    RawFiles = [];
    FileFormat = [];
end
% Force list of files to be a cell array
if ~isempty(RawFiles) && ischar(RawFiles)
    RawFiles = {RawFiles};
end
% Some verifications
if ~isempty(RawFiles) && isempty(FileFormat)
    error('If you pass the filenames in input, you must define also the FileFormat argument.');
end
% Get Protocol information
ProtocolInfo = bst_get('ProtocolInfo');
OutputDataFile = [];


%% ===== SELECT DATA FILE =====
% If file to load was not defined : open a dialog box to select it
if isempty(RawFiles) 
    % Get default import directory and formats
    LastUsedDirs = bst_get('LastUsedDirs');
    DefaultFormats = bst_get('DefaultFormats');
    % Get MRI file
    [RawFiles, FileFormat, FileFilter] = java_getfile( 'open', ...
        'Open raw EEG/MEG recordings...', ...  % Window title
        LastUsedDirs.ImportData, ...           % Last used directory
        'multiple', 'files_and_dirs', ...      % Selection mode
        bst_get('FileFilters', 'raw'), ...     % List of available file formats
        DefaultFormats.DataIn);
    % If no file was selected: exit
    if isempty(RawFiles)
        return
    end
    % Save default import directory
    LastUsedDirs.ImportData = bst_fileparts(RawFiles{1});
    bst_set('LastUsedDirs', LastUsedDirs);
    % Save default import format
    DefaultFormats.DataIn = FileFormat;
    bst_set('DefaultFormats',  DefaultFormats);
    % Process the selected directories :
    %    1) If they are .ds/ directory with .meg4 and .res4 files : keep them as "files to open"
    %    2) Else : add all the data files they contains (subdirectories included)
    RawFiles = io_expand_filenames(FileFilter, RawFiles);
    if isempty(RawFiles)
        bst_error(['No data ' FileFormat ' file in the selected directories.'], 'Open raw EEG/MEG recordings...', 0);
        return
    end
end


%% ===== IMPORT =====
iOutputStudy = [];
% Loop on the files to import
for iFile = 1:length(RawFiles)
    % ===== OPENING FILE =====
    bst_progress('start', 'Open raw EEG/MEG recordings', 'Reading file header...');
    % Open file
    [sFile, ChannelMat, errMsg] = in_fopen(RawFiles{iFile}, FileFormat, ImportOptions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Update ChannelMat %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('ChannelMat.mat');
    megidx = find(strcmp({ChannelMat.Channel.Type}, 'MGG')==1); %Does nothing at the moment
    %sensorNames149 = {ChannelMat.Channel.Name};
    
    for i = 1:length(megidx)
        sensorName = ChannelMat.Channel(megidx(i)).Name;
        
        ChannelMat.Channel(megidx(i)).Loc = ChannelMat.Channel(megidx(i)).Loc;

        %ChannelMat.Channel(megidx(i)).Loc = ChannelMat.Channel(ismember(sensorNames149,sensorName)).Loc;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(sFile)
        bst_progress('stop');
        return;
    end
    % Yokogawa non-registered warning
    if ~isempty(errMsg) && ImportOptions.DisplayMessages
        java_dialog('warning', errMsg, 'Open raw EEG/MEG recordings');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Convert to Continues Data %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Messages = [];
%     % Check if loaded file is epoched
%         if isempty(sFile.epochs) || (length(sFile.epochs) == 1) || strcmpi(sFile.format, 'CTF-CONTINUOUS')
%             Messages = 'Only the files that contain two epochs or more can be converted to continuous files.';
%             sFile = [];
%             return;
%         end
%         % Process each epoch
%         for i = 1:length(sFile.epochs)
%             % Rebuild absolute sample indices for this file
%             nSamples = sFile.epochs(i).samples(2) - sFile.epochs(i).samples(1) + 1;
%             if (i == 1)
%                 epochSmp = [0, nSamples-1];
%             else
%                 epochSmp = [epochSmp; epochSmp(i-1,2) + [1, nSamples]];
%             end
%         end
%         % Update events
%         for iEvt = 1:length(sFile.events)
%             % Detect events that appear a the last sample of a trial and at the first one of the next one (only for simple events)
%             if (size(sFile.events(iEvt).times,1) == 1)
%                 % Get the length of the epoch in samples for each event occurrence
%                 smpEpoch = [sFile.epochs(sFile.events(iEvt).epochs).samples];
%                 % Detect if the occurrence is at the first sample or in the last 5 samples
%                 isLast  = (smpEpoch(2:2:end) - sFile.events(iEvt).samples < 5);
%                 isFirst = (sFile.events(iEvt).samples - smpEpoch(1:2:end) == 0);
%                 % Detect the markers that are doubled: last sample of epoch #i and first of epoch #i+1 
%                 iDouble = find(isFirst(2:end) & isLast(1:end-1)) + 1;
%             else
%                 iDouble = [];
%             end
%             % Process each occurrence of each independent separately
%             for iOcc = 1:size(sFile.events(iEvt).samples,2)
%                 % Identify the epoch in which this event is occurring
%                 iEpoch = sFile.events(iEvt).epochs(iOcc);
%                 adjustSmp = epochSmp(iEpoch,1) - sFile.epochs(iEpoch).samples(1);
%                 % Re-refence this event occurrence starting from the beginning of the continuous file (sample 0)
%                 sFile.events(iEvt).samples(:,iOcc) = sFile.events(iEvt).samples(:,iOcc) + adjustSmp;
%             end
%             % Update times and epoch indice
%             sFile.events(iEvt).times  = sFile.events(iEvt).samples ./ sFile.prop.sfreq;
%             sFile.events(iEvt).epochs = ones(size(sFile.events(iEvt).epochs));
%             % Remove those doubled markers (remove the first sample of epoch #i+1)
%             if ~isempty(iDouble)
%                 % Get the times to remove
%                 tRemoved = sFile.events(iEvt).times(1,iDouble);
%                 % Remove the events occurrences
%                 sFile.events(iEvt).times(:,iDouble)   = [];
%                 sFile.events(iEvt).samples(:,iDouble) = [];
%                 sFile.events(iEvt).epochs(:,iDouble)  = [];
%                 if ~isempty(sFile.events(iEvt).reactTimes)
%                     sFile.events(iEvt).reactTimes(:,iDouble) = [];
%                 end
%                 % Display message
%                 Messages = [Messages, 10, 'Removed ' num2str(length(iDouble)) ' x "' sFile.events(iEvt).label, '": ', sprintf('%1.3fs ', tRemoved)];
%             end
%         end
%         % Display message
%         if ~isempty(Messages)
%             Messages = ['Errors detected in the events of the .ds file (duplicate markers): ' Messages];
%         end
%         
%         % Remove epochs 
%         sFile.epochs = [];
%         % Update the other fields
%         sFile.prop.samples = [epochSmp(1,1), epochSmp(end,2)];
%         sFile.prop.times   = sFile.prop.samples ./ sFile.prop.sfreq;
%         sFile.format = 'CTF-CONTINUOUS';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % ===== OUTPUT STUDY =====
    % Get short filename
    [tmp, fBase] = bst_fileparts(RawFiles{iFile});
    % Build output condition name
    if isfield(sFile, 'condition') && ~isempty(sFile.condition)
        ConditionName = ['@raw' sFile.condition];
    else
        ConditionName = ['@raw' fBase];
    end
    % Output subject
    if isempty(iSubject)
        % Get default subject
        SubjectName = 'NewSubject';
        [sSubject, iSubject] = bst_get('Subject', SubjectName, 1);
        % If subject does not exist yet: create it
        if isempty(sSubject)
            [sSubject, iSubject] = db_add_subject(SubjectName);
        end
        % If subject cannot be created
        if isempty(sSubject)
            bst_error(['Could not create subject "' SubjectName '"'], 'Open raw EEG/MEG recordings');
            return;
        end
    else
        % Get specified subject
        sSubject = bst_get('Subject', iSubject, 1);
    end
    % Do not allow automatic registration with head points when using the default anatomy
    if (sSubject.UseDefaultAnat)
        ImportOptions.ChannelAlign = 0;
    end

    % If condition already exists
    [sExistStudy, iExistStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, file_standardize(ConditionName, 1)));
    if ~isempty(sExistStudy) && ~isempty(sExistStudy.Data)
        % Need to check if the raw file is the same or they are two files with the same name in different folders
        % Get the raw data files
        iRaw = find(strcmpi({sExistStudy.Data.DataType}, 'raw'));
        if ~isempty(iRaw)
            % Load data description
            DataFile = sExistStudy.Data(iRaw).FileName;
            DataMat = in_bst_data(DataFile);
            % If same filenames: cannot link it again in the database
            LinkFile = DataMat.F.filename;
            minLength = min(length(LinkFile), length(RawFiles{iFile}));
            if file_compare(LinkFile(1:minLength), RawFiles{iFile}(1:minLength))
                %bst_error('This file is already available in the explorer.', 'Open raw EEG/MEG recordings', 0);
                panel_protocols('SelectNode', [], 'rawdata', iExistStudy, iRaw );
                OutputDataFile = DataFile;
                continue;
            % Else: Create a condition with a different name
            else
                % Add a numeric tag at the end of the condition name
                curPath = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sExistStudy.FileName));
                curPath = file_unique(curPath);
                [tmp__, ConditionName] = bst_fileparts(curPath, 1);
            end
        end
    end
    % Create output condition
    iOutputStudy = db_add_condition(sSubject.Name, ConditionName);
    if isempty(iOutputStudy)
        error('Folder could not be created : "%s/%s".', bst_fileparts(sSubject.FileName), ConditionName);
    end
    % Get output study
    sOutputStudy = bst_get('Study', iOutputStudy);
    % Get the study in which the channel file has to be saved
    [sChannel, iChannelStudy] = bst_get('ChannelForStudy', iOutputStudy);
    
    % ===== CREATE DEFAULT CHANNEL FILE =====
    if isempty(ChannelMat)
        ChannelMat = db_template('channelmat');
        ChannelMat.Comment = [sFile.device ' channels'];
        ChannelMat.Channel = repmat(db_template('channeldesc'), [1, length(sFile.channelflag)]);
        % For each channel
        for i = 1:length(ChannelMat.Channel)
            if (length(ChannelMat.Channel) > 99)
                ChannelMat.Channel(i).Name = sprintf('E%03d', i);
            else
                ChannelMat.Channel(i).Name = sprintf('E%02d', i);
            end
            ChannelMat.Channel(i).Type    = 'EEG';
            ChannelMat.Channel(i).Loc     = [0; 0; 0];
            ChannelMat.Channel(i).Orient  = [];
            ChannelMat.Channel(i).Weight  = 1;
            ChannelMat.Channel(i).Comment = [];
        end
        % Add channel file to database
        db_set_channel(iChannelStudy, ChannelMat, 0, 0);

    % ===== SAVE LOADED CHANNEL FILE =====
    else
        % Add history field to channel structure
        ChannelMat = bst_history('add', ChannelMat, 'import', ['Link to file: ' RawFiles{iFile} ' (Format: ' FileFormat ')']);
        % Remove fiducials only from polhemus and ascii files
        isRemoveFid = ismember(FileFormat, {'MEGDRAW', 'POLHEMUS', 'ASCII_XYZ', 'ASCII_NXYZ', 'ASCII_XYZN', 'ASCII_NXY', 'ASCII_XY', 'ASCII_NTP', 'ASCII_TP'});
        % Perform the NAS/LPA/RPA registration for some specific file formats
        isAlign = ismember(FileFormat, {'NIRS-BRS'});
        % Detect auxiliary EEG channels
        ChannelMat = channel_detect_type(ChannelMat, isAlign, isRemoveFid);
        % Do not align data coming from Brainstorm exported files (already aligned)
        if strcmpi(FileFormat, 'BST-BIN')
            ImportOptions.ChannelAlign = 0;
        end
        % Add channel file to database
        [ChannelFile, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Modality] = db_set_channel(iChannelStudy, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign);
        % Display the registration if this was skipped in db_set_channel
        if (ImportOptions.ChannelAlign == 0) && (ImportOptions.DisplayMessages) && ~isempty(Modality)
            bst_memory('UnloadAll', 'Forced');
            channel_align_manual(ChannelFile, Modality, 0);
        end
    end
    
    % ===== SAVE DATA FILE =====
    % Build output filename
    OutputDataFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sOutputStudy.FileName), ['data_0raw_' fBase '.mat']);
    % Build output structure
    DataMat = db_template('DataMat');
    DataMat.F           = sFile;
    DataMat.Comment     = 'Link to raw file';
    DataMat.ChannelFlag = sFile.channelflag;
    DataMat.Time        = sFile.prop.times;
    DataMat.DataType    = 'raw';
    DataMat.Device      = sFile.device;
    % Compumedics: add start time to the file comment
    if strcmpi(sFile.format, 'EEG-COMPUMEDICS-PFS') && isfield(sFile.header, 'rda_startstr') && ~isempty(sFile.header.rda_startstr)
        DataMat.Comment = [DataMat.Comment ' [' sFile.header.rda_startstr ']'];
    end
    % Add history field
    DataMat = bst_history('add', DataMat, 'import', ['Link to raw file: ' RawFiles{iFile}]);
    % Save file on hard drive
    bst_save(OutputDataFile, DataMat, 'v6');
    % Add file to database
    sOutputStudy = db_add_data(iOutputStudy, OutputDataFile, DataMat);

    % ===== UPDATE DATABASE =====
    % Update links
    db_links('Study', iOutputStudy);
    % Refresh both data node and channel node
    iUpdateStudies = unique([iOutputStudy, iChannelStudy]);
    panel_protocols('UpdateNode', 'Study', iUpdateStudies);
end

% If something was imported
if ~isempty(iOutputStudy)
    % Select the data study node
    panel_protocols('SelectStudyNode', iOutputStudy);
    % Save database
    db_save();
end

bst_progress('stop');
end 