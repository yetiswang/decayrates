
directories = dir();

dirFlags = [directories.isdir];
subFolders = directories(dirFlags);


pcles = length(subFolders) - 2; 

for i = 3 : length(subFolders)
    cd(subFolders(i).name)
    load decayrates
    [PCR] = PhotonCountRate( decayrates );
    close all
    clear decayrates
    cd .. 
end 
