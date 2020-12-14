function out_data_egi_raw( DataMat, OutputFile)
% OUT_DATA_EGI_RAW: Exports a MEG/EEG recordings file in EGI NetStation RAW format
%
% USAGE:  out_data_egi_raw( DataMat, OutputFile );
%
% INPUT: 
%    - DataMat    : full path to Brainstorm file to export
%    - OutputFile : full path to output file (with '.raw' extension)

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2014 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2008

% Define precision
precision = 'real*4';
precisionCode = 4;

% Open file
fid = fopen(OutputFile, 'wb', 'b');
if (fid < -1)
    error(['Could not open file "' OutputFile '".']);
end

%% ===== WRITE HEADER =====
% Get current time
c = clock();
% Get sampling rate
samplingRate = (1/(DataMat.Time(2) - DataMat.Time(1)));
% RAW header
fwrite(fid, precisionCode, 'integer*4');  % versionNumber :  (2:integer*4, 4:real*4, 6:real*8)
fwrite(fid, c(1),          'integer*2'); % Year
fwrite(fid, c(2),          'integer*2'); % Month
fwrite(fid, c(3),          'integer*2'); % Day
fwrite(fid, 0,             'integer*2'); % Hour
fwrite(fid, 0,             'integer*2'); % Minute
fwrite(fid, 0,             'integer*2'); % Second
fwrite(fid, 0,             'integer*4'); % Milisecond
% fwrite(fid, samplingRate,      'double'); % samplingRate % Not needed in
% 2019
fwrite(fid, samplingRate,      'integer*2'); % samplingRate
fwrite(fid, size(DataMat.F,1), 'integer*2'); % numChans
fwrite(fid, 1,                 'integer*2'); % boardGain
fwrite(fid, 0,                 'integer*2'); % numConvBits
fwrite(fid, 0,                 'integer*2'); % ampRange
fwrite(fid, size(DataMat.F,2), 'integer*4'); % numSamples
fwrite(fid, 1,                 'integer*2'); % numEvents
fwrite(fid, 'tim0',            'uchar');     % (time = 0)

%% ===== WRITE SAMPLES =====
% Write sample / sample (all the channels for each time frame)
for iTime = 1:size(DataMat.F, 2)
    % Values for all channels
    fwrite(fid, DataMat.F(:,iTime) * 1e6, precision);
    % Value for event 'tim0'
    if (abs(DataMat.Time(iTime)) < 1e-5) 
        val = 1;
    else
        val = 0;
    end
    fwrite(fid, val, precision);
end

% Close file
fclose(fid);



