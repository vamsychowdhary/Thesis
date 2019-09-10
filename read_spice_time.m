function [time,amp]=read_spice_time(fname)

% This function reads an input file in .txt format containing the data
% format as that of the output file form LTSpice. Only reads transient
% domain files.

% Open the file
fid=fopen(fname);

% Read header
linecounter=0;
data = [];
myStop=[];

while isempty(myStop)
    tline = fgetl(fid);
    linecounter = linecounter+1;
    myStop =  str2num(tline(1)); % waits for a number
end
% Header is finished when first character in line is a number 

% Read data
while ~isempty(myStop)
    linedata=sscanf(tline,'%f %f',[1,2]);
    data= [data; linedata];
    tline = fgetl(fid);
    linecounter = linecounter+1;
    if feof(fid)
        myStop=[];  % EOF
    else
        myStop =  feof(fid) & str2num(tline(1)); % waits for a number
    end
end

% Close file
fclose(fid);

% Output data
time=data(:,1);
amp=data(:,2);
return;