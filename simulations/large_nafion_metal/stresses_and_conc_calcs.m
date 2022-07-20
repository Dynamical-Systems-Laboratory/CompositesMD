%
% Script for computation of average RDFs and CNs
%

% Directories to consider
dir_names = dir('seed_*');
dir_names = {dir_names.name};
ndirs = length(dir_names);
nbins = 50;

% All data
% Stresses and contributions
all_s_xx = zeros(ndirs, nbins);
all_s_yy = zeros(ndirs, nbins);
all_s_zz = zeros(ndirs, nbins);
all_s_xy = zeros(ndirs, nbins);
all_s_yz = zeros(ndirs, nbins);
all_s_xz = zeros(ndirs, nbins);
all_s_ke_xx = zeros(ndirs, nbins);
all_s_virial_xx = zeros(ndirs, nbins);
all_s_ke_yy = zeros(ndirs, nbins);
all_s_virial_yy = zeros(ndirs, nbins);
all_s_ke_zz = zeros(ndirs, nbins);
all_s_virial_zz = zeros(ndirs, nbins);
all_s_kinetic_energy = zeros(ndirs, nbins);
all_s_potential_energy = zeros(ndirs, nbins);

for i=1:ndirs
   fname = dir_names(i);
   path = sprintf('%s',fname{:});
   % Move to path and execute the python script for averaging
   cd(sprintf('%s/post_processing/', path));

   % Collect the data
   temp = load('stresses_and_eng_out_xx.txt');
   all_s_xx(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_yy.txt');
   all_s_yy(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_zz.txt');
   all_s_zz(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_xy.txt');
   all_s_xy(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_yz.txt');
   all_s_yz(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_xz.txt');
   all_s_xz(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_ke_xx.txt');
   all_s_ke(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_virial_xx.txt');
   all_s_virial(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_kinetic_energy.txt');
   all_s_kinetic_energy(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_potential_energy.txt');
   all_s_potential_energy(i,:) = mean(temp);
   % yy
   temp = load('stresses_and_eng_out_ke_yy.txt');
   all_s_ke_yy(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_virial_yy.txt');
   all_s_virial_yy(i,:) = mean(temp);
   % zz
   temp = load('stresses_and_eng_out_ke_zz.txt');
   all_s_ke_zz(i,:) = mean(temp);
   temp = load('stresses_and_eng_out_virial_zz.txt');
   all_s_virial_zz(i,:) = mean(temp);
        
   cd '../../'
end
% Save data
save('stress_conc_data');


