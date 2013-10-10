close all
clear all

p = genpath('/data/scripts/DCE/Loveless2011_201204MULTI');%/home/thomasn/scripts/ReferenceRegionRatios');%/home/data/scripts/DCE/ReferenceRegion2/ReferenceRegionYankRatio');
path(p,path);

CPU = 80;
CHECK = 0;
r2filter = 0;
timestamp = 45;
xdata = [];

NAMES = { 'tn14_24feb12', 'tn15_24feb12', 'tn16_24feb12', 'tn01_24feb12', 'tn10_24feb12', 'tn04_24feb12'};%,'tn02_20dec11', 'tn04_20dec11', 'tn11_20dec11', 'tn12_20dec11'};%, 'tn12_20dec11'};%, 'rj14_07mar12'};%,'tn01_24feb12',};%tn02_29dec11','tn01_29dec11', 'tn24_20dec11', 'tn22_20dec11'};%, 'rj17_15nov11', 'rj08_15nov11', 'rj04_15nov11','rj18_15nov11', 'rj15_15nov11','rj03_15nov11','rj04_15nov11','tn05_08sep11', 'tn09_08sep11', 'tn13_08sep11', 'tn14_08sep11', 'tn15_08sep11', 'tn19_08sep11'}%, 'rj04_11aug11', 'rj10_11aug11', 'rj03_11aug11', ,};

%DATES(1).isdir = 1; DATES(1).name = '20120113';
%DATES(2).isdir = 1; DATES(2).name = '20111222';
%DATES(1).isdir = 1; DATES(1).name = '20111223';

place = '/data/studies/PETMRI';
for i = 1:numel(NAMES)
    D = dir(fullfile(place, NAMES{i}));
    
         
%                 if(~isempty(DATES))
%                     D = DATES;
%                 end
%     
    for j = 1:numel(D)
        
        if(D(j).isdir)
            if(isempty(strfind(D(j).name, '.')))
                
                % Date
            
                
                G = dir(fullfile(place, NAMES{i}, D(j).name));
           
                    
                %%FIX!!!
                for qq = 1:numel(G)
                    
                    if(G(qq).isdir)
                        
                        if(isempty(strfind(G(qq).name, '.')))
                            
                            
                            % type
                            E = dir(fullfile(place, NAMES{i}, D(j).name, G(qq).name));
                            
                            for k = 1:numel(E)
                                
                                
                                if(E(k).isdir)
                                    if(isempty(strfind(E(k).name, '.')))
                                        if(~isempty(strfind(E(k).name, 'dynamic')))
                                            
                                            F = dir(fullfile(place, NAMES{i}, D(j).name, G(qq).name,E(k).name));
                                            
                                            
                                            
                                            for l = 1:numel(F)
                                                
                                                if(~(F(l).isdir))
                                                    if(~isempty(strfind(F(l).name, 'fitted_R1info_RRratios.mat')))
                                                        
                                                        %fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        
                                                        % D_fit_RRFXLvoxelsC(fullfile(place, NAMES{i}, D(j).name, G(qq).name,E(k).name, F(l).name), CPU, CHECK, r2filter)
                                                     %   C_fit_AIFnovp_voxels(fullfile(place, NAMES{i}, D(j).name, G(qq).name,E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                        %disp('Moo')
                                                    elseif(~isempty(strfind(F(l).name, '01_R1info.mat')))
                                                        disp('Cool')
                                                        fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        
                                                    %  D_fit_AIFwithvp_voxels_COH(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                   xdata{end+1}.out = D_fit_AUC_RR_voxels_COH(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter, timestamp);
                                                    xdata{end}.name = NAMES{i};
                                                    xdata{end}.study= E(k).name;
                                                    elseif(~isempty(strfind(F(l).name, '01_AIF_FXR_ROI.mat')))
                                                      %  disp('Fool')
                                                       % fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        
                                                      % D_fit_FXR_voxels_COH(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter);
                                                    elseif(~isempty(strfind(F(l).name, 'AIF_iRGDfull_ROI.mat')))
                                                        %                                                         disp('iRGD')
                                                        %                                                         fullfile(place, NAMES{i}, D(j).name, E(k).name, F(l).name)
                                                        %
                                                        %  C_fit_iRGD_voxels(fullfile(place, NAMES{i}, D(j).name, G(qq).name, E(k).name, F(l).name), CPU, CHECK, r2filter);
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
