function varargout = process_ra_beta_import_sara_data( varargin )
%R. Avci May 2018
%macro_methodcall;
eval(macro_method)
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Read Data';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'SARA';
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
    %OutputFiles{1} = ra_test_import_raw(FileName, FileFormat, iSubject, ImportOptions);
    OutputFiles{1} = import_raw(FileName, FileFormat, iSubject, ImportOptions);
end