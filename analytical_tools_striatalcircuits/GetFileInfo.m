function [file, neurons, events] = GetFileInfo(file);
% GetFileInfo.m  Laubach Lab: Eyal Kimchi, Kumar Narayanan 2006
% This code loads the events and neurons from a .nex (NeuroExplorer) file
% Usage:  [file, neurons, events] = GetFileInfo(file), where file is a
% struct, and the file name is specified as file.name.  neurons and events
% are lists of chars of neuron and event names.  For use with nex_ts, a
% plexon utility for loading variables out of nex files. 
% Example: 
% file.name = 'myfile.nex'; 
% [file, neurons, events] = GetFileInfo(file)
% note that the file struct will be overwritten with file details

% Counter to open the file
file.start = 0;
% assume 1 hour to start
file.stop = 3600;
file.moredata = [];

fid = fopen(file.name, 'r');    % Opens the file
if(fid == -1)
	fprintf('Can not open .nex file.\n');
   return
end
magic = fread(fid, 1, 'int32');         %Reads the file
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');         %Frequency of digitization
file.start = fread(fid, 1, 'int32')/freq;
file.stop = fread(fid, 1, 'int32')/freq;
nvar = fread(fid, 1, 'int32');          % Number of variables in the header
fseek(fid, 260, 'cof');
names = zeros(1, 64);                   %Skips past the header
for i=1:nvar
	types(i) = fread(fid, 1, 'int32');
	var_version = fread(fid, 1, 'int32');
	names(i, :) = fread(fid, [1 64], 'char');
	dummy = fread(fid, 128+8, 'char');
end
names = char(names);
fclose(fid);

var_names = cellstr(names);             % Converts the variable names to a string
neurons = strmatch('SPK', names);       % Gets the names of all the neurons
neurons = var_names(neurons);
unique_length = 7; % for plexon signals, the character value that is unique
% [neurons, file.num_neurons] = GetSignalNames(neurons, unique_length);   %A function that finds unique neurons
events = strmatch('EV', names);
events = var_names(events);
file.num_events = length(events);


function [unique_names, num_unique] = GetSignalNames(names, unique_length);
% GetFileInfo.m  Laubach Lab: Eyal Kimchi, Kumar Narayanan 2006
% This code finds which neurons are unique
% Usage: [unique_names, num_unique] = GetSignalNames(names, unique_length);
% Where names is a list of names, and unique_length is the character that is unique. 
% Example: names = ['fun1', 'fun2']; unique_length=4;
% [unique_names, num_unique] = GetSignalNames(names, unique_length);


if isempty(names) == 1
    unique_names = []; num_unique = 0; 
    return;
else
	for i = 1:length(names)
        crop_names{i} = names{i}(1:unique_length);
	end
	
	unique_names = unique(crop_names);
	num_unique = length(unique_names);
end