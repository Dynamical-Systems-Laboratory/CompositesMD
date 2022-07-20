% Computes average Ep of the system during the last 3 ns NVT
% Average is over all seeds

% Directories to consider
dir_names = dir('seed_*');
dir_names = {dir_names.name};
ndirs = length(dir_names);

% Collect the data
all_seeds = zeros(ndirs, 1);
for i=1:ndirs
    fname = dir_names(i);
    path = sprintf('%s/',fname{:});

    % Move to path and load data
    cd(sprintf('%s', path));
    fprintf("Processing %s\n", path)
    
    data = load('thermo_data.txt');

    all_seeds(i) = mean(data(size(data,1)-3000:end,5));

    cd '../'
end

mean(all_seeds)
std(all_seeds)
