% Script for plotting all the stress components (need to do it manually)

clear; close all;
load('stress_conc_data')

% Optionally plot only some seeds
seeds_to_plot = [1, 10];

% Number of seeds
nseeds = 9;

% Colors (depend on the atoms, preferences)
clr_back = [ 128, 249, 91; 253, 155, 152; 175, 175, 175; 137, 182, 249]/255;
clr_front = [0, 85, 0; 165, 0, 0; 0, 0, 0; 0, 0, 127]/255;

% Box and bin dimensions
box = [4.7486326432291293e+01 2.1549589316776195e+02;
4.8789955212853897e-02 6.7402386444758392e+01;
1.6272271435415320e-01 6.7288453685675051e+01];
nbins = 50;

dx=box(1,2)-box(1,1);
L=dx/10;
x=0:L/nbins:L;
xc=x(2:end)-L/nbins/2;

kplt = 1;
for i=1:nseeds
      
    % For all of them
    subplot(3, 3, i);
    % xx
    h = plot(xc, all_s_xx(i,:), 'Color', clr_front(1,:), 'LineWidth', 2);
    % yy
%     h = plot(xc, all_s_yy(i,:), 'Color', clr_front(2,:), 'LineWidth', 2);
    % zz
%     h = plot(xc, all_s_zz(i,:), 'Color', clr_front(3,:), 'LineWidth', 2);
    
    grid on
    set(gca,'FontSize',14,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
        'on');

    xlim([0,20])

% For x
    ylim([-6e3, 6e3])

    yticks([-6e3, -3e3, 0, 3e3, 6e3])
    yticklabels({'-6', '-3', '0', '3', '6'})


% For y
%     ylim([-2e4, 2e4])
% 
%     yticks([-2e4, -1e4, 0, 1e4, 2e4])
%     yticklabels({'-2', '-1', '0', '1', '2'})

% For z
%     ylim([-2e4, 2e4])
% 
%     yticks([-2e4, -1e4, 0, 1e4, 2e4])
%     yticklabels({'-2', '-1', '0', '1', '2'})

 % Create ylabel
    ylabel('','Interpreter','latex');

    % Create xlabel
    xlabel('$x$ direction, [nm]','Interpreter','latex');

end

