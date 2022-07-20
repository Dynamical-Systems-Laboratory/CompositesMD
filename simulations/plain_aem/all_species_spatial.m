% Script for plotting all the number densities

clear; close all;

% Colors (depend on the atoms, preferences)
clr = [0,0,1; 0,0,1; 0,0,1; 0,0,1; 0,0,1; 1,0,0; 0,0,1; 0,0,0; ...
        0.5,0.5,0.5; 0.5,0.5,0.5];

% Directories to consider
dir_names = dir('seed_*');
dir_names = {dir_names.name};
ndirs = length(dir_names);

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
        sp=ff(ix0(2)+1:ixf-1)     % return the substring up to 2nd underscore

        data = load(ff);
        temp = mean(data, 2);
        all_data(inm,:) = temp;
%         plot(temp/max(max(temp)), 'Color', clr(inm,:), 'LineWidth', 2)
        plot(temp, 'Color', clr(inm,:), 'LineWidth', 2)
        species{end+1} = sp;
        hold on
    end
    legend(species)
    cd '../../'
end
