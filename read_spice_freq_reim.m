function [freq,mag,phase]=read_spice_freq_reim(fname)

% This function reads an input file in .txt format containing the data
% format as that of the output file form LTSpice. Please remember to select
% cartesian (re,im) format before saving the file.

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
    linedata = replace(tline,","," ");
    linedata=sscanf(linedata,'%f %f %f',[1 3]);
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
freq=data(:,1);
mag=20*log10(sqrt(data(:,2).^2+data(:,3).^2));
phase=atan2(data(:,3),data(:,2));
return;