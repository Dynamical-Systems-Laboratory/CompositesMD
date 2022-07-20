%% Post-process RDF averages
clear, close all
load('stress_conc_data')

% Box and bin dimensions in A (from seed 1)
box = [3.7896407571804986e+01 2.3277265642818821e+02;
3.1399400417377166e-01 7.2918558095847203e+01;
3.5368210001783495e-01 7.2878869999947568e+01];
% Bin width is 3.8975 A
nbins = 50;

dx=box(1,2)-box(1,1);
L=dx/10;
x=0:L/nbins:L;
xc=x(2:end)-L/nbins/2;

% Volume of one bin, m3
V_bin = dx/nbins*(box(2,2)-box(2,1))*(box(3,2)-box(3,1))*1e-30;
% Temperature, K
T = 300;
% Universal gas constant, atm*m3/(mol K)
R = 8.2057e-5;
% Avogadro's number
N_avg = 6.0221409e+23;

plot_sc_error_area(xc, all_s_xx, 1, [137, 182, 249]/255, [0, 0, 127]/255, '$\sigma_{xx}$, [atm]')
plot_sc_error_area(xc, all_s_yy, 2, [253, 155, 152]/255, [165, 0, 0]/255, '$\sigma_{yy}$, [atm]')
plot_sc_error_area(xc, all_s_zz, 3, [128, 249, 91]/255, [0/255, 85/255, 0/255], '$\sigma_{zz}$, [atm]')
% plot_sc_error_area(xc, all_s_virial_xx, 'virial xx', 4, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
% plot_sc_error_area(xc, all_s_ke_xx, 'kinetic xx', 5, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
% 
%plot_sc_error_area(xc, all_s_kinetic_energy, 'Kinetic energy', 6, [186/255, 149/255, 193/255],[134/255, 5/255, 159/255], 'Eenergy, [Kcal/mol]')
plot_sc_error_area(xc, all_s_potential_energy, 7, [253, 155, 152]/255, [165, 0, 0]/255, 'Energy, [Kcal/mol]')

function plot_sc_error_area(x, y, i, clrB, clrF, ylab)
    
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
    
    fill([x;flipud(x)],[y-dy;flipud(y+dy)], clrB*0.8,'linestyle','none');
    plot(x,y, 's-', 'LineWidth',2, 'Color', clrF*0.5, 'MarkerSize',  10)

    hold off
    
    % Create ylabel
    ylabel(ylab,'Interpreter','latex');

    % Create xlabel
    xlabel('x direction, [nm]','Interpreter','latex');

    % Create title
%     title(plot_title,'Interpreter','latex');

    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim(axes1,[0 5]);
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
        'on');
   
    set(gcf, 'Position', [417   289   790   582])
    
end