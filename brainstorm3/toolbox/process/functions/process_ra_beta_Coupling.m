function varargout = process_ra_beta_Coupling(varargin)
% process_test_ra plots fHR
% @=======================================================================
%macro_methodcall;
eval(macro_method)
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()%#ok<DEFNU>
        % Description the process
    sProcess.Comment     = 'Analysis-Coupling';
    sProcess.FileTag     = 'PSD';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'SARA';
    sProcess.Index       = 1301;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};%, 'results', 'timefreq', 'matrix'};;
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)%#ok<DEFNU>
    Comment = ['Running ' sProcess.Comment ' Analysis'];
end

%% ===== RUN =====
function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    sOutput = {};
end