% Script for caculating the ratio of interface and bulk thickenss

clear; close all;

% Colors (depend on the atoms, preferences)
clr_back = [ 128, 249, 91; 253, 155, 152; 175, 175, 175; 137, 182, 249]/255;
clr_front = [0, 85, 0; 165, 0, 0; 0, 0, 0; 0, 0, 127]/255;

% Plotting groups (lists of indices in all_data that correspond to atoms in
% that group)
polymer = [1:2, 4, 5, 7]; 
ion = [3];
platinum = [6];
water = [8, 9];

% Box and bin dimensions
box = [3.7896407571804986e+01 2.3277265642818821e+02;
3.1399400417377166e-01 7.2918558095847203e+01;
3.5368210001783495e-01 7.2878869999947568e+01];
nbins = 50;

dx=box(1,2)-box(1,1);
L=dx/10;
Lbin = L/nbins;
x=0:L/nbins:L;
xc=x(2:end)-L/nbins/2;
    
% Directories to consider
dir_names = dir('seed_*');
dir_names = {dir_names.name};
ndirs = length(dir_names);

% All data
all_interface = zeros(ndirs, 1);
all_bulk = zeros(ndirs, 1);
for i=1:ndirs
    fname = dir_names(i);
    path = sprintf('%s/post_processing/',fname{:});

    % Move to path and load data
    cd(sprintf('%s', path));
    fprintf("Processing %s\n", path)

    file_names = dir('number_density_*');
    species = {};
    for inm = 1:length(file_names)
        ff = file_names(inm).name;
        ix0=strfind(ff,'_');  % get the underscore locations
        ixf=strfind(ff,'.');  % get the . location
        sp=ff(ix0(2)+1:ixf-1);     % return the substring up to 2nd underscore

        temp = load(ff);
        all_data(i,:,inm) = mean(temp, 2);
    end
    
    ps = sum(all_data(i,:,polymer),3);
    is = sum(all_data(i,:,ion),3);
    pts = sum(all_data(i,:,platinum),3);
    ws = sum(all_data(i,:,water),3);
    
    % Locate the bulk (Pt = 0)
    ib_all = find(pts == 0);
    ib0 = min(ib_all);
    ibF = max(ib_all);
    
    % Bulk length
    all_bulk(i) = Lbin*length(ib_all);
    
    % Interface width (total)
    % All indices of bins that make up the interface - right and left
    int_all = find((pts ~= 0) & ((is ~= 0) | (ps ~= 0) | (ws ~= 0)));
    % Technically this is not needed - just multiply the size of previous 
%     part = find(int_all(1:end-1) ~= int_all(2:end) - 1);
    all_interface(i) = Lbin*length(int_all);

    cd '../../'
end

mean(all_interface./all_bulk)
std(all_interface./all_bulk)
