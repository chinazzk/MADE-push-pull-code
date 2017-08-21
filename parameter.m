%*************************************
% Coded by Hayden Fischer, Spring 2014
%
% Copies parameters from Ucode output
% and places them in SLIM tcl files
%*************************************

cd ../..

% Read parameter.txt
fid = fopen('parameter_ADE.txt');
par = textscan(fid, '%f %*[^\n]');

phi     = par{1}(1);
alpha_l = par{1}(2);
alpha_t = par{1}(3);

cd SLIM_Facies/ADE


%*****************
% Edit inject file
%*****************

% Read porosity line
fid = fopen('madesite_slim_inject.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set phi ',8);
    counter = counter + 1;
end
line_phi = counter;

% Read longitudinal dispersivity line
fid = fopen('madesite_slim_inject.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_l',11);
    counter = counter + 1;
end
line_alpha_l = counter;

% Read transverse dispersivity line
fid = fopen('madesite_slim_inject.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_t',11);
    counter = counter + 1;
end
line_alpha_t = counter;

% Read entire inject file
fid = fopen('madesite_slim_inject.tcl');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

% Rewrite inject file
A{line_phi}     = char(strcat({'set alpha_l              '},num2str(alpha_l)));
A{line_alpha_l} = char(strcat({'set alpha_t              '},num2str(alpha_t)));
A{line_alpha_t} = char(strcat({'set phi                  '},num2str(phi)));

fid = fopen('madesite_slim_inject.tcl', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end


%***************
% Edit wait file
%***************

% Read porosity line
fid = fopen('madesite_slim_wait.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set phi ',8);
    counter = counter + 1;
end
line_phi = counter;

% Read longitudinal dispersivity line
fid = fopen('madesite_slim_wait.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_l',11);
    counter = counter + 1;
end
line_alpha_l = counter;

% Read transverse dispersivity line
fid = fopen('madesite_slim_wait.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_t',11);
    counter = counter + 1;
end
line_alpha_t = counter;

% Read entire wait file
fid = fopen('madesite_slim_wait.tcl');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

% Rewrite wait file
A{line_phi}     = char(strcat({'set alpha_l              '},num2str(alpha_l)));
A{line_alpha_l} = char(strcat({'set alpha_t              '},num2str(alpha_t)));
A{line_alpha_t} = char(strcat({'set phi                  '},num2str(phi)));

fid = fopen('madesite_slim_wait.tcl', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end


%******************
% Edit extract file
%******************

% Read porosity line
fid = fopen('madesite_slim_extract.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set phi ',8);
    counter = counter + 1;
end
line_phi = counter;

% Read longitudinal dispersivity line
fid = fopen('madesite_slim_extract.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_l',11);
    counter = counter + 1;
end
line_alpha_l = counter;

% Read transverse dispersivity line
fid = fopen('madesite_slim_extract.tcl');
trip    = 0;
counter = 0;
while trip == 0
    tline = fgets(fid);
    trip = strncmp(tline,'set alpha_t',11);
    counter = counter + 1;
end
line_alpha_t = counter;

% Read entire extract file
fid = fopen('madesite_slim_extract.tcl');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

% Rewrite extract file
A{line_phi}     = char(strcat({'set alpha_l              '},num2str(alpha_l)));
A{line_alpha_l} = char(strcat({'set alpha_t              '},num2str(alpha_t)));
A{line_alpha_t} = char(strcat({'set phi                  '},num2str(phi)));

fid = fopen('madesite_slim_extract.tcl', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

exit