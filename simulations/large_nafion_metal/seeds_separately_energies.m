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
    h = plot(xc, all_s_potential_energy(i,:), 'Color', clr_front(4,:), 'LineWidth', 2);

 set(gca,'FontSize',14,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
        'on');

    grid on
    xlim([0,20])
    ylim([-7e4, 0])

    yticks([-7e4, -3.5e4, 0])
    yticklabels({'-7', '-3.5', '0'})

    % Create ylabel
    ylabel('','Interpreter','latex');

    % Create xlabel
    xlabel('$x$ direction, [nm]','Interpreter','latex');
end

