function varargout = process_bl_test_create_raw_data_v2( varargin )

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Import SQUID data';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'MGG';
    sProcess.Index       = 999;
    sProcess.isSeparator = 0;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = [];
    sProcess.nMinFiles   = 0;
        
    % Option: Subject name
    sProcess.options.subjectname.Comment = 'Please click run to import data';
    sProcess.options.subjectname.Type    = 'label';
    sProcess.options.subjectname.Value   = 'Please click run to import data';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

%% Use GUI to find filepaths
c = import_MGG_txt_app; % Output from GUI are squid array and squid data filepaths (squidarray_filepath and squiddata_filepath)
waitfor(c);
clear c
load filepaths
%% Squid positions

% Run if you have squid position data otherwise comment this section out
[arraySquid, unmodifiedSensorNames] = readData(squidarray_filepath); % Contains sensor names + locations

% Create channelmat file needed for brainstorm import function
ChannelMat=CreateSensorMat(arraySquid.Name, arraySquid.Type, arraySquid.posData, arraySquid.oriData);

save ChannelMat.mat ChannelMat

%% Data input and formatting

% Import data 
A = readmatrix(squiddata_filepath);
sensorName = (A(1, 2:end));
sensorName = sensorName(~isnan(sensorName));
Time = A(2:end,1);
mgg_data = A(2:end, 2:length(sensorName) + 1);
Calibration = calib_setting;
SenIdsUsed = [33 3 35 9 15 25 23 26 30 36 5 7 16 17 12 14 18 22 19];
OutPath = './';
RawFile = 'BS_RawFile';


if length(sensorName) < length(unmodifiedSensorNames)
            x = ismember(unmodifiedSensorNames, string(sensorName));
            mgg_all = zeros(length(Time), length(unmodifiedSensorNames));
            mgg_all(:, x) = mgg_data;
            mgg_data = mgg_all;
end

if Calibration
    
    CalibFac = readmatrix(calib_filepath);
    [~, Calid] = ismember(SenIdsUsed, CalibFac(:,1));
    CalX = CalibFac(Calid,2);
    
    % First find the indices of mgg_data that correspond to SenIdsUsed and then
    % arrange the data in the order of SenIdsUsed in a new array
    for i = 1:length(SenIdsUsed)
        
        % Iterate through SenIdsUsed one by one and find index of match in the
        % sensor names of the header file
        [logical, index] = ismember((SenIdsUsed(i)), sensorName);
        if logical > 0
            % Arrange mgg_data in the order of SenIdsUsed in a new array
            mgg_data_in_order_of_senID(:, i) = mgg_data(:, index);
        end
    end
    
    % Fill new array mgg_data_cal with data from SenIdsUsed that has been
    % calibrated with a calibration factor
    for i = 1:length(SenIdsUsed)
        
        % Iterate through SenIdsUsed one by one and find index of match in ALL
        % sensor names
        [logical, index] = ismember(string(SenIdsUsed(i)), unmodifiedSensorNames);
        if logical > 0
            % If a match is found, populate an array that will contain ALL the
            % mgg_data (including calibrated data). Iterate through
            % mgg_data_in_order to calibrate data gathered by SenIdsUsed.
            
            % It is important to populate an array with ALL 41 mgg_data
            % otherwise it will not be compatible with the 41 channel.mat data
            mgg_data_cal(:, index) = mgg_data_in_order_of_senID(:, i) .* repmat(CalX(i),1,size(mgg_data_in_order_of_senID(:,i),1))';
        end
    end
    
    % Fill in remaining mgg_data in mgg_data_cal pertaining to sensors that
    % haven't been calibrated.
    for i = 1:length(sensorName)
        
        % Iterate through sensor names from the header but only look at sensors
        % that have not been used in SenIdsUsed.
        if ~ismember(sensorName(i), SenIdsUsed)
            
            % Find index of sensor in header of text file that matches the current
            % sensor name. If a match is found, populate the column with
            % original mgg_data
            [logical, index] = ismember(string(sensorName(i)), unmodifiedSensorNames);
            if logical > 0
                mgg_data_cal(:, index) = mgg_data(:, i); % i represent the indices of mgg_data that need to be scaled by calibration factor
            end
        end
    end
    
    RawFile = 'BS_RawFile_Calibrated';
    mgg_data = mgg_data_cal;
   
end

%% Create RAW brainstorm file

OutFile=strcat(OutPath,RawFile, '.raw');
Comment = 'MGG';
% Set the correct sampling frequency
sf = 1/(Time(2)-Time(1)); %
% Assing the data in the correct format
datamat.data=mgg_data;
datamat.sf=sf;
datamat.Comment=Comment;
[datamatbs]=CreatedbBrainstorm(datamat,1:size(mgg_data,2));
out_data_egi_raw(datamatbs, OutFile);


end

%% Helper function for finding SQUID locations and names 
function [dataSquid, unmodifiedSensorNames]=readData(filePath)

% Import data 
A = importdata(filePath,'\t',1); % Assumes .txt file has only 1 header and is tab delimited
xPos = A.data(:,4);
yPos = A.data(:,5);
zPos = A.data(:,6) + 5;

xOri = A.data(:,7);
yOri = A.data(:,8);
zOri = A.data(:,9);

unmodifiedSensorNames = cellstr(string(num2cell(A.data(:,1)')));

sensorNames = cellstr(string(num2cell(A.data(:,1)')));
sensorNames_REF = sensorNames(1:12);
sensorNames_nonREF = sensorNames(13:end);
sensorType = A.textdata(2:end, 3);

centre = cellstr(string(num2cell([24, 33, 34])));
inner_circle = cellstr(string(num2cell([3, 9, 15, 23, 25, 35])));
outer_circle = cellstr(string(num2cell([4, 5, 7, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 26, 27, 30, 36, 37])));

one = cellstr(string(num2cell([22, 23])));
two = cellstr(string(num2cell([19])));
three = cellstr(string(num2cell([3, 20, 26, 27])));
four = cellstr(string(num2cell([30])));
five = cellstr(string(num2cell([35, 36])));
six = cellstr(string(num2cell([4, 5, 10])));
seven = cellstr(string(num2cell([7, 9])));
eight = cellstr(string(num2cell([16])));
nine = cellstr(string(num2cell([11, 13, 15, 17])));
ten = cellstr(string(num2cell([12])));
eleven = cellstr(string(num2cell([14, 25])));

xMagsensor = cellstr(string(num2cell([10, 13, 20, 21, 34])));
yMagsensor = cellstr(string(num2cell([4, 11, 24, 27, 37])));

for i = 1:length(sensorNames_REF)
    sensorNames_REF{i} = ['REF_sensor ' sensorNames_REF{i}];
end 

for i = 1:length(sensorNames_nonREF)
    if ismember(sensorNames_nonREF{i},xMagsensor)
        coord = '_x';
    elseif ismember(sensorNames_nonREF{i},yMagsensor)
        coord = '_y';
    else
        coord = '';
    end
    
    if ismember(sensorNames_nonREF{i},centre)
        ringLetter = 'A';
        clockNumber = '0';
    end 
    if ismember(sensorNames_nonREF{i},inner_circle) 
        ringLetter = 'B';
        if ismember(sensorNames_nonREF{i},one)
            clockNumber = '1';
        elseif ismember(sensorNames_nonREF{i},two)
            clockNumber = '2';
        elseif ismember(sensorNames_nonREF{i},three)
            clockNumber = '3';
        elseif ismember(sensorNames_nonREF{i},four)
            clockNumber = '4';
        elseif ismember(sensorNames_nonREF{i},five)
            clockNumber = '5';
        elseif ismember(sensorNames_nonREF{i},six)
            clockNumber = '6';
        elseif ismember(sensorNames_nonREF{i},seven)
            clockNumber = '7';
        elseif ismember(sensorNames_nonREF{i},eight)
            clockNumber = '8';
        elseif ismember(sensorNames_nonREF{i},nine)
            clockNumber = '9';
        elseif ismember(sensorNames_nonREF{i},ten)
            clockNumber = '10';
        elseif ismember(sensorNames_nonREF{i},eleven)
            clockNumber = '11';
        else
            clockNumber = '12';
        end 
    end
    if ismember(sensorNames_nonREF{i},outer_circle)
        ringLetter = 'C';
        if ismember(sensorNames_nonREF{i},one)
            clockNumber = '1';
        elseif ismember(sensorNames_nonREF{i},two)
            clockNumber = '2';
        elseif ismember(sensorNames_nonREF{i},three)
            clockNumber = '3';
        elseif ismember(sensorNames_nonREF{i},four)
            clockNumber = '4';
        elseif ismember(sensorNames_nonREF{i},five)
            clockNumber = '5';
        elseif ismember(sensorNames_nonREF{i},six)
            clockNumber = '6';
        elseif ismember(sensorNames_nonREF{i},seven)
            clockNumber = '7';
        elseif ismember(sensorNames_nonREF{i},eight)
            clockNumber = '8';
        elseif ismember(sensorNames_nonREF{i},nine)
            clockNumber = '9';
        elseif ismember(sensorNames_nonREF{i},ten)
            clockNumber = '10';
        elseif ismember(sensorNames_nonREF{i},eleven)
            clockNumber = '11';
        else
            clockNumber = '12';
        end 
    end
    sensorNames_nonREF{i} = [ringLetter clockNumber coord];
end 

sensorNames = [sensorNames_REF, sensorNames_nonREF]';
dataSquid.Name = sensorNames;
dataSquid.Type = sensorType;
dataSquid.posData = [xPos, yPos, zPos];
dataSquid.oriData = [xOri, yOri, zOri];


%% Uncomment to check sensor locations with SQUID map

% scatter3(xPos,yPos,zPos);
% for i = 1:length(sensorName)
%     a = sensorName(i);
%     b = num2str(a); % will need to change this if I'm already using a string obviously
%     c = cellstr(b);
%     dx = 1;
%     dy = 1;
%     dz = 1;% displacement so the text does not overlay the data points
%     text(xPos(i)+dx, yPos(i)+dy, zPos(i) + dz, c);
% end

end

%%
function [ChannelMat]=CreateSensorMat(sensorName, typeData, posData, oriData)

ChannelMat = db_template('channelmat');
Channel = ChannelMat.Channel;

for m=1:size(posData,1)
    Loc = zeros(3,8);
    Loc(:,1)=posData(m,:)';
    Loc(:,2)=posData(m,:)';
    Channel(1,m).Loc=Loc;
    Channel(1,m).Orient=oriData(m, :)';
    Channel(1,m).Comment='Vanderbilt Magnetometer';
    Channel(1,m).Group='';
    Channel(1,m).Weight=[0.25;0.25;0.25;0.25;-0.25;-0.25;-0.25;-0.25];
    Channel(1,m).Type= char(typeData(m));
    
    Channel(1,m).Name = char(sensorName(m));
end
ChannelMat.Channel = Channel;
end 