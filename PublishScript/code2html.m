files = dir('hw_*.m');
for id = 1 : length(files)
    publish(files(id).name);
    close all;
end