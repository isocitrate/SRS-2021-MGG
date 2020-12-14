
function [datamat]=CreatedbBrainstorm(datamat,iGoodChannels)


% TimeAll=(0:size(dch,1)-1)/sf;
% sizeBlock=180*sf; 
% nBlocks = ceil(length(TimeAll) / sizeBlock);
% i=1;%wich block of trial
% ind = min([(i-1)*sizeBlock + 1, i*sizeBlock], length(TimeAll)); 
% data=dch(ind(1):ind(2),:);
% Time = TimeAll(ind(1):ind(2));
% figure, plot(Time,data');
% xlabel('Time (s)');ylabel('fT')

data=datamat.data;
sf=datamat.sf;
filename=datamat.Comment;
if isfield(datamat,'Time')
  Time=datamat.Time;
else
  Time=(0:size(data,1)-1)/sf;
end

nSensors = size(data,2);
datamat = db_template('datamat');
datamat.Comment = filename;
datamat.Device='CTF';
datamat.Time = Time;
datamat.ChannelFlag = repmat(-1,nSensors,1);
datamat.ChannelFlag(iGoodChannels)=1;
datamat.F = zeros(nSensors, length(datamat.Time));
datamat.F = data';


end