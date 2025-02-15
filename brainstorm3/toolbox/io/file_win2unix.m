function FileName = file_win2unix( FileName )
% FILE_WIN2UNIX: Transform a Matlab filename to universal filename (replaces '\' with '/').

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
% Authors: Francois Tadel, 2009

% Nothing in input
if isempty(FileName)
    return
end
% Replace backslashes
FileName = strrep(FileName, '\', '/');
% Remove first slash
%if (FileName(1) == '/') && ~file_exist(FileName)
%   FileName(1) = [];

if (isequal(FileName(1), '/')) && (~file_exist(FileName))
    FileName(1) = [];
    
end



