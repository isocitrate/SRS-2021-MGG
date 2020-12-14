function varargout = process_ra_beta_HRv(varargin)
% process_test_ra plots fHR
% @=======================================================================
%macro_methodcall;
eval(macro_method)
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()%#ok<DEFNU>
        % Description the process
    sProcess.Comment     = 'Analysis-HRv';
    sProcess.FileTag     = '| HRv';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'SARA';
    sProcess.Index       = 1302;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};%, 'results', 'timefreq', 'matrix'};;
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Definition of the options
    sProcess.options.mMCG.Comment = 'mHRv';
    sProcess.options.mMCG.Type    = 'checkbox';
    sProcess.options.mMCG.Value   = 0;
    sProcess.options.fMCG.Comment = 'fHRv';
    sProcess.options.fMCG.Type    = 'checkbox';
    sProcess.options.fMCG.Value   = 0;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)%#ok<DEFNU>
    ismMCG = sProcess.options.mMCG.Value;
    isfMCG = sProcess.options.fMCG.Value;
    if ismMCG && ~isfMCG
        Comment = ['Running ' sProcess.Comment ' Analysis for mMCG'];
    elseif isfMCG && ~ismMCG
        Comment = ['Running ' sProcess.Comment ' Analysis for fMCG'];
    elseif ~ismMCG && ~isfMCG
        Comment = 'Select mMCG or fMCG';
    elseif ismMCG && isfMCG
        Comment = 'Select either mMCG or fMCG';
    else
        Comment = '';
    end
end

%% ===== RUN =====
function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    clc;
    clear sOutput;
    sOutput = {};
    ismMCG = sProcess.options.mMCG.Value;
    isfMCG = sProcess.options.fMCG.Value;
    % Get Data from BS
    DataMat = in_bst_data(sInputs.FileName);
    % Get Events from BS or compute
    events = DataMat.F.events;
    sf = length(DataMat.Time)/DataMat.Time(end);
    lbls = {events.label}';
    if ismMCG
        event_id = find(strcmp(lbls,'mMCG')==1);
    elseif isfMCG
        event_id = find(strcmp(lbls,'fMCG')==1);
    end
    iR = events(event_id).samples;
    rran=get_rr_intervals(iR,sf);%[samples x 1]    
    SDs_poincare=get_poincareHRV(rran);
    % Calculate HRV measures
    [mRR, ~, RMSsd, ApEn, pNN20, pNN15, pNN10] = stats(rran);
    % HR
    [HR, time]=get_HR(iR,sf);
    isplot = 1;
    if isplot
        %FigHandle1=plotPoincare(rran,SDs_poincare,sInputs.FileName); 
        FigHandle2=plotHRate(iR,sf,strcat(sInputs.FileName),2);
    end
    
    % ===== SAVE THE RESULTS =====
    % Get the output study (pick the one from the first file)
    iStudy = sInputs(1).iStudy;
    % Create a new data file structure
    chflags = (-1)*abs(DataMat.F.channelflag);
    chflags(1) = 1;
    DataMat2 = db_template('datamat');
    DataMat2.F           = repmat(HR,length(chflags),1);
    if ismMCG
        DataMat2.Comment     = 'mHRv';
    elseif isfMCG
        DataMat2.Comment     = 'fHRv';
    end
    DataMat2.ChannelFlag = chflags;   % List of good/bad channels (1=good, -1=bad)
    DataMat2.Time        = time;
    DataMat2.DataType    = 'recordings';
    DataMat2.nAvg        = 1;
    % Create a default output filename
    sOutput{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'data_custom');
    % Save on disk
    save(sOutput{1}, '-struct', 'DataMat2');
    % Register in database
    db_add_data(iStudy, sOutput{1}, DataMat2);
    
    % Display results on MATLAB
    disp('Heart Rate Variability Metrics')
    disp(['mean HR:' num2str(mean(HR)) ' bpm']);
    disp(['mean RR:' num2str(1000*mRR) ' ms']);
    disp(['RMSsd:' num2str(1000*RMSsd) ' ms']);
    disp(['ApEn:' num2str(ApEn)]);
    disp(['pNN20:' num2str(pNN20)])
    disp(['pNN15:' num2str(pNN15)])
    disp(['pNN10:' num2str(pNN10)])
    disp(['poincare-SD1:' num2str(SDs_poincare.SD1) ' ms']);
    disp(['poincare-SD1:' num2str(SDs_poincare.SD2) ' ms']);
end