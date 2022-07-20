% Script for caculating the stress oscillations in the bulk and the average
% stress in the bulk

clear; close all;

% Plotting groups (lists of indices in all_data that correspond to atoms in
% that group)
polymer = [1:5, 7]; 
ion = [6];
platinum = [8];
water = [9, 10];
    
% Box and bin dimensions
box = [4.7486326432291293e+01 2.1549589316776195e+02;
4.8789955212853897e-02 6.7402386444758392e+01;
1.6272271435415320e-01 6.7288453685675051e+01];
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

% Load the stress data
main_data = load('stress_conc_data.mat');

% All data (xx, yy, and zz components)
all_min = zeros(ndirs, 3);
all_max = zeros(ndirs, 3);
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
    
    % Minimum and maximum stress values
    all_min(i,1) = min(main_data.all_s_xx(i, ib0:ibF));
    all_min(i,2) = min(main_data.all_s_yy(i, ib0:ibF));
    all_min(i,3) = min(main_data.all_s_zz(i, ib0:ibF));
    
    all_max(i,1) = max(main_data.all_s_xx(i, ib0:ibF));
    all_max(i,2) = max(main_data.all_s_yy(i, ib0:ibF));
    all_max(i,3) = max(main_data.all_s_zz(i, ib0:ibF));

    cd '../../'
end

disp('Average max-min')
fprintf("%f, %f, %f\n" , mean(all_max(:,1)-all_min(:,1)), ...
    mean(all_max(:,2)-all_min(:,2)), ...
    mean(all_max(:,3)-all_min(:,3)))

disp('Std max-min')
fprintf("%f, %f, %f\n" , std(all_max(:,1)-all_min(:,1)), ...
    std(all_max(:,2)-all_min(:,2)), ...
    std(all_max(:,3)-all_min(:,3)))

disp('Max max-min')
fprintf("%f, %f, %f\n" , max(all_max(:,1)-all_min(:,1)), ...
    max(all_max(:,2)-all_min(:,2)),...
    max(all_max(:,3)-all_min(:,3)))

disp('Min max-min')
fprintf("%f, %f, %f\n" , min(all_max(:,1)-all_min(:,1)), ...
    min(all_max(:,2)-all_min(:,2)), ...
    min(all_max(:,3)-all_min(:,3)))


