% Script for plotting all the number densities

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
x=0:L/nbins:L;
xc=x(2:end)-L/nbins/2;
    
% Directories to consider
dir_names = dir('seed_*');
dir_names = {dir_names.name};
ndirs = length(dir_names);

% All data
all_elastic = zeros(ndirs, 1);
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
    
    ps = sum(all_data(i,:,polymer),3);
    is = sum(all_data(i,:,ion),3);
    pts = sum(all_data(i,:,platinum),3);
    ws = sum(all_data(i,:,water),3);
    
    % For all of them
    subplot(3, 3, i); 
    h = plot(xc, ps/max(ps),  xc, is/max(is), xc, pts/max(pts), xc, ws/max(ws), 'LineWidth', 2);
          xlim([0,20])
        ylim([0,1])
    % For a subset
%     subplot(2, 1, kplt);
%     if ismember(i, seeds_to_plot)
%         h = plot(xc, ps/max(ps),  xc, is/max(is), xc, pts/max(pts), xc, ws/max(ws), 'LineWidth', 2);
%         kplt = kplt + 1;
%     end
%     

  %  axes1 = axes('Parent',gcf);
   set(gca,'FontSize',14,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
        'on');

    for ic = 1:length(clr_front)
        set(h(ic), 'color', clr_front(ic,:));
    end
    grid on
%     axis tight
%     set(gca, 'fontname', 'arial')
    set(gca, 'fontsize', 14)
   
        % Create ylabel
    ylabel('','Interpreter','latex');

    % Create xlabel
    xlabel('$x$ direction, [nm]','Interpreter','latex');
    
%     legend(species)
    cd '../../'
end

