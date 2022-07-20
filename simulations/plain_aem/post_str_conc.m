%% Post-process RDF averages
clear, close all
load('stress_conc_data')

% Box and bin dimensions
box = [9.3164746561401799e+00 1.2217463514386839e+02;
3.1786516451474611e-01 3.3407723035481752e+01;
-3.5382023519630934e-02 3.3760970223506611e+01];
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

plot_sc_error_area(xc, all_s_xx, 'xx', 1, [186/255, 149/255, 193/255],[134/255, 5/255, 159/255], 'Stress, [atm]')
plot_sc_error_area(xc, all_s_yy, 'yy', 2, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
plot_sc_error_area(xc, all_s_zz, 'zz', 3, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
% plot_sc_error_area(xc, all_s_virial_xx, 'virial xx', 4, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
% plot_sc_error_area(xc, all_s_ke_xx, 'kinetic xx', 5, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Stress, [atm]')
% 
% plot_sc_error_area(xc, all_s_kinetic_energy, 'Kinetic energy', 6, [186/255, 149/255, 193/255],[134/255, 5/255, 159/255], 'Eenergy, [Kcal/mol]')
% plot_sc_error_area(xc, all_s_potential_energy, 'Potential energy', 7, [180/255, 209/255, 223/255],[17/255, 122/255, 175/255], 'Energy, [Kcal/mol]')

function plot_sc_error_area(x, y, plot_title, i, clrB, clrF, ylab)
    
    % Create figure
    figure1 = figure(i);

    % Create axes
    axes1 = axes('Parent',figure1);
    
    % Initial
    x = x';
    
    hold on
    
    % Final
%     dy = std(y)';
%     y = mean(y)';
    
%     fill([x;flipud(x)],[y-dy;flipud(y+dy)], clrB*0.8,'linestyle','none');
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
   
end