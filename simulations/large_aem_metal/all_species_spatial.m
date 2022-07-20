% Script for plotting all the number densities

clear; close all;

% Colors (depend on the atoms, preferences)
clr_back = [ 128, 249, 91; 253, 155, 152; 175, 175, 175; 137, 182, 249]/255;
clr_front = [0, 85, 0; 165, 0, 0; 0, 0, 0; 0, 0, 127]/255;

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
x=0:L/nbins:L;
xc=x(2:end)-L/nbins/2;

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
        sp=ff(ix0(2)+1:ixf-1);     % return the substring up to 2nd underscore

        temp = load(ff);
        all_data(i,:,inm) = mean(temp, 2);
%         plot(temp/max(max(temp)), 'Color', clr(inm,:), 'LineWidth', 2)
%         plot(temp, 'Color', clr(inm,:), 'LineWidth', 2)
%         species{end+1} = sp;
%         hold on
    end
%     legend(species)
    cd '../../'
end

% Plot the groups
temp =  sum(all_data(:,:,polymer),3);
plot_sc_error_area(xc,temp./max(temp')', 'Distribution of atoms', 1, clr_back(1,:), clr_front(1,:), 'Backbone', 'atoms')

temp = sum(all_data(:,:,ion),3);
add_to_plot(xc, temp./max(temp')', 'Distribution of atoms', 1, clr_back(2,:), clr_front(2,:), 'Counterion', 'atoms')

temp = sum(all_data(:,:,platinum),3);
add_to_plot(xc, temp./max(temp')', 'Distribution of atoms', 1, clr_back(3,:), clr_front(3,:), 'Platinum', 'atoms')

temp = sum(all_data(:,:,water),3);
add_to_plot(xc, temp./max(temp')', 'Distribution of atoms', 1, clr_back(4,:), clr_front(4,:), 'Water', 'atoms')

function plot_sc_error_area(x, y, plot_title, i, clrB, clrF, ylab, figname)
    
    % Create figure
    figure1 = figure(i);

    % Create axes
    axes1 = axes('Parent',figure1);
    
    % Initial
    x = x';
    
    hold on
    
    % Final
    dy = std(y)';
    y = mean(y)';
    
    fill([x;flipud(x)],[y-dy;flipud(y+dy)], clrB*0.8,'linestyle','none', 'FaceAlpha',.3,'EdgeAlpha',.3);
    plot(x,y, 's-', 'LineWidth',2, 'Color', clrF*0.5, 'MarkerSize',  10)

    hold off
    
    % Create ylabel
    ylabel(ylab,'Interpreter','latex');

    % Create xlabel
    xlabel('x direction, [nm]','Interpreter','latex');

    % Create title
    title(plot_title,'Interpreter','latex');

    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim(axes1,[0 5]);
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
        'on');
   
    set(gcf, 'Position', [417   289   790   582])
    
    savefig(figname)
    
end

function add_to_plot(x, y, plot_title, i, clrB, clrF, ylab, figname)
    
    close all

    figure1 = openfig(figname);


    % Create axes
%     axes1 = axes('Parent',figure1);
    
    % Initial
    x = x';
    
    hold on
    
    % Final
    dy = std(y)';
    y = mean(y)';
    
    fill([x;flipud(x)],[y-dy;flipud(y+dy)], clrB*0.8,'linestyle','none', 'FaceAlpha',.3,'EdgeAlpha',.3);
    plot(x,y, 's-', 'LineWidth',2, 'Color', clrF*0.5, 'MarkerSize',  10)

    hold off
    
    % Create ylabel
%     ylabel(ylab,'Interpreter','latex');
% 
%     % Create xlabel
%     xlabel('x direction, [nm]','Interpreter','latex');
% 
%     % Create title
%     title(plot_title,'Interpreter','latex');
% 
%     % Uncomment the following line to preserve the Y-limits of the axes
%     % ylim(axes1,[0 5]);
%     box(axes1,'on');
%     % Set the remaining axes properties
%     set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
%         'on');
   
%     set(gcf, 'Position', [417   289   790   582])
    
    savefig(figname)
end