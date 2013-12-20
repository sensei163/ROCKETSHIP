close all
clear all

p = genpath('/home/thomasn/scripts');
path(p,path);

CPU = 60;
CHECK = 0;
r2filter = 0;

NAMES = { 'rj03_15nov11', 'rj04_11aug11', 'rj10_11aug11', 'rj03_11aug11'};
place = '/home/thomasn/data';
for i = 1:numel(NAMES)
    D = dir(fullfile(place, NAMES{i}));
    
    for j = 1:numel(D)
        
        if(D(j).isdir)
            if(isempty(strfind(D(j).name, '.')))
                
                
                
                G = dir(fullfile(place, NAMES{i}, D(j).name));
                
                for qq = 1:numel(G)
                    
                    if(G(qq).isdir)
                        
                        if(isempty(strfind(G(qq).name, '.')))
                            E = dir(fullfile(place, NAMES{i}, D(j).name, G(qq).name));
                            
                            for k = 1:numel(E)
                                
                                
                                if(E(k).isdir)
                                    if(isempty(strfind(E(k).name, '.')))
                                        if(~isempty(strfind(E(k).name, 'dynamic')))
                                            
                                            F = dir(fullfile(place, NAMES{i}, D(j).name, G(qq).name,E(k).name));
                                            
                                            
                                            for l = 1:numel(F)
                                                
                                                if(~(F(l).isdir))
                                                    if(~isempty(strfind(F(l).name, 'newAIF_no_vpFIT_ROI.mat')))
                                                        
                                                        %fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        %C_fit_AIFnovp_voxels(fullfile(place, NAMES{i}, D(j).name, G(qq).name,E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                        %disp('Moo')
                                                    elseif(~isempty(strfind(F(l).name, 'newAIF_with_vpFIT_ROI.mat')))
                                                        disp('Cool')
                                                         fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        
                                                       D_fit_AIFwithvp_voxels_COH(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                    elseif(~isempty(strfind(F(l).name, 'newAIF_FXR_ROI.mat')))
                                                       disp('Fool')
                                                        fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        
                                                      D_fit_FXR_voxels_COH(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                    else
                                                        
                                                        
                                                        
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
end
