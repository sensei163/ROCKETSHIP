function [meanAIF, meanSignal] = AIF_auto_cluster( concentration_array, image_array, time_vector, TR, species ) 
%  This function, MEAN_AIF, is used to generate an AIF automatically via K-means clustering. 
% 
% concentration_array is simply the image data in matrix form. 
% 
% MEAN_AIF is outputted AIF. 
% }

%%%%%%%%%%%%%%%%%%%%%% %FUNCTION BODY BEGINS HERE%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checking the size of the image array: Is it a slice or a stack?
image_dims = ndims(concentration_array); 
if image_dims == 4   % a stack
[dimx, dimy, dimz,dimt] = size(concentration_array); 
end 
if image_dims == 3   % a slice 
[dimx,dimy,dimt] = size(concentration_array); 
dimz = 1; 
end 

% Doing some housekeeping: Defining Variables 
ctcs = zeros(dimx*dimy*dimz,dimt); 
time_2_peak = zeros(dimx*dimy*dimz,1);       % Time to peak of concentration time course (CTC)
peak_height = zeros(dimx*dimy*dimz,1);       % Peak height, signal value, arb. units. 
bolus_time = zeros(dimx*dimy*dimz,1);        % Time of contrast agent injection. 
ctc_integral = zeros(dimx*dimy*dimz,1);      % Time integral of the CTC.  
fm_integral = zeros(dimx*dimy*dimz,1);       % First moment of time integral. 
fm_peak_ratio = zeros(dimx*dimy*dimz,1);       % First moment of time integral. 


% Eventually we want to correlate the indices of the above arrays with
% x,y,and z coordinates. This is done below: 

index_array = zeros(dimx*dimy*dimz,3); 

% Cleaning the image array: replacing bad values with workable ones. 

concentration_array(concentration_array<0) = 0;        % concentrations should not be negative. 
concentration_array(isnan(concentration_array)) = 0;   % NaNs just screw things up. Make them 0. 
concentration_array(concentration_array>32767) = 0;    % Infs screw things up too. This gets rid of them 

%% Select voxels and compute stuff about them: 

% IMPORTANT: If the data is bad in a certain aspect it will be obliterated.
% So too will all the correspending quantities.  This will take place in
% the loop below.  

jj = 1; 
for k = 1 : dimz 
    for j = 1 : dimy 
        for i = 1 : dimx
            
            if image_dims == 3
                
                ctc = squeeze(concentration_array(i,j,:)); 
                 
            else
                ctc = squeeze(concentration_array(i,j,k,:));  
            end
            
            if max(ctc) > 0 %correponding end at line 164 
                
            if ctc(1) > 0 
                ctc = ctc - ctc(1);
                ctc(find(ctc<0)) = 0; 
            end 
           
            % This is an array of time courses. 
            ctcs(jj,:) = ctc; 
            
            % Filling our image array. 
            index_array(jj,:) = [i,j,k];  
            
           
           
            % Now, Perform Computations (on each voxel) 
            
            %Peak height and time-to-peak can be outputted by max
            %function:
            [peak_height(jj),time_2_peak(jj)] = max(ctc); 
            
           %Taking the derivative; a useful quantity. 
           
           ctc_deriv = diff(ctc); 
          
           % Now we will compute the ctc time integral. We will do this via
           % trapezoid rule. 
           
           %What time would you like to integrate according to: 
           times = 0 : length(ctc)-1;
           times = times'; 
           
           % These integrals are 'supposed' to be conducted on
           % concentration time courses, Thus we are going to set the
           % minimum value of each the curves to zero and go from there. 
           
           %ctc = ctc - min(ctc);  
           
           %Now integrate: 
           ctc_integral(jj) = trapz(ctc); 
           
           % Onto the ctc first moment. We integrate different points: 
           % concenration (signal value) multiples by the time, t; 
           
           ctc_FM = ctc .* times; 
           
           fm_integral(jj) = trapz(ctc_FM); 
           
           % Now the FM_Peak_ratio:
           
           fm_peak_ratio(jj) = peak_height(jj)/ fm_integral(jj); 
           
           % Now we inspect the values and remove the bad ones:  
           
           %Checking for smoothness: Using the second derivative. 
           spikiness = abs(mean(ctc_deriv)); 
           smoothness = diff(ctc_deriv); 
           smoothness = abs(mean(smoothness)); % One value and we want it absolute to make things easier. 
           
           % Removing voxels that have ackward ctc. If they look nothing
           % like a ctc, remove them! 
           if strcmp(species, 'human')
               if max(ctc) <= 1.5 || smoothness >= 0.05 || max(ctc)>100 || spikiness > 5 || time_2_peak(jj) > 5
                   ctcs(jj,:) = [];
                   time_2_peak(jj) = [];       
                   peak_height(jj) = [];        
                   bolus_time(jj) = [];         
                   ctc_integral(jj) = [];      
                   fm_integral(jj) = [];  
                   index_array(jj,:) = []; 
                   fm_peak_ratio(jj) = []; 

                   jj = jj -1; 

               end
           elseif strcmp(species, 'mouse') 
                 if max(ctc) <= 1.5 || smoothness >= 2 || max(ctc)>200 || spikiness > 10 || time_2_peak(jj) > 8
                   ctcs(jj,:) = [];
                   time_2_peak(jj) = [];       
                   peak_height(jj) = [];        
                   bolus_time(jj) = [];         
                   ctc_integral(jj) = [];      
                   fm_integral(jj) = [];  
                   index_array(jj,:) = []; 
                   fm_peak_ratio(jj) = []; 

                   jj = jj -1; 

                 end
           end
               
           
            else
               ctcs(jj,:) = [];
               time_2_peak(jj) = [];       
               peak_height(jj) = [];        
               bolus_time(jj) = [];         
               ctc_integral(jj) = [];      
               fm_integral(jj) = [];  
               index_array(jj,:) = []; 
               fm_peak_ratio(jj) = []; 
               
               jj = jj -1; 
            end
               
                
         jj = jj + 1; 
        end
    end
end
                
%% Cluster Analysis: 

% At this point we need to arrange our row-vectors into an array that can
% be sent to the KMEANS function.  

% While we calculated a bunch of quantities, lets try sorting with only
% two: the FM integral and the peak height.  

clustering_variables = cat(2,peak_height,fm_integral,fm_peak_ratio,time_2_peak); %,time_2_peak,ctc_integral); 
clustering_variables(clustering_variables == inf) = 0; 
[clusters,~,~,~] = kmeans(clustering_variables,5,'Replicates', 50, 'EmptyAction', 'drop','Distance', 'cityblock'); 

%% Define the AIF:

% The AIF should be associated with the cluster that has the lowest FM
% integral values and highest peak height values.   

% We are going to plot them and ask the used which cluster they want to use
% as the AIF: 

good_aif = 0; 

while good_aif ~= 1 
close all;

% Lets plot the curves from the different clusters: 
figure; 

%time vector for the plots below 
if ndims(concentration_array) == 4  
    tx = 0 : length(concentration_array(1,1,1,:))-1; 
    tx = tx .* TR; 
else 
    tx = 0 : length(concentration_array(1,1,:))-1; 
    tx = tx .* TR; 
end
    

subplot(2,4,1); 
plot(tx,ctcs(clusters==1,:),'r'); 
xlabel('Time'); 
ylabel('Concentration'); 
hold off; 
subplot(2,4,2); 
plot(tx,ctcs(clusters==2,:),'g'); 
xlabel('Time'); 
ylabel('Concentration'); 
hold off; 
subplot(2,4,3); 
plot(tx,ctcs(clusters==3,:),'b'); 
xlabel('Time'); 
ylabel('Concentration'); 
hold off; 
subplot(2,4,4); 
plot(tx,ctcs(clusters==4,:),'m'); 
xlabel('Time'); 
ylabel('Concentration'); 
hold off; 
subplot(2,4,5); 
plot(tx,ctcs(clusters==5,:),'c'); 
xlabel('Time'); 
ylabel('Concentration'); 
hold off; 


%% Plotting the clusters onto the original image to show locations: 

% find the elements in each clusuter: 
    %Cluster 1:
    cluster_1_indices = find(clusters==1);  
    cluster_2_indices = find(clusters==2); 
    cluster_3_indices = find(clusters==3); 
    cluster_4_indices = find(clusters==4); 
    cluster_5_indices = find(clusters==5); 
    
    % making the array (in 3 and 4D)
    if image_dims == 3
    rgb_cluster_array = zeros(dimy,dimx,3); 
        for i = 1 : length(clustering_variables(:,1)) 
            if find(cluster_1_indices == i) 
                %Make these voxels red 
                rgb_cluster_array(index_array(i,2),index_array(i,1),1) = 1; 
            elseif find(cluster_2_indices == i) 
                rgb_cluster_array(index_array(i,2),index_array(i,1),2) = 1; 
            elseif find(cluster_3_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),3) = 1;
            elseif find(cluster_4_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),1) = 1;
                rgb_cluster_array(index_array(i,2),index_array(i,1),3) = 1;
            elseif find(cluster_5_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),2) = 1;
                rgb_cluster_array(index_array(i,2),index_array(i,1),3) = 1;
            end
        end
    figure;
    image(rgb_cluster_array); 
    axis image; 
    hold off; 
    axis off; 
  
    else  % Now the 4D stack case...
    rgb_cluster_array = zeros(dimy,dimx,3,dimz); 
    
        for i = 1 : length(clustering_variables(:,1)) 
            if find(cluster_1_indices == i) 
                %Make these voxels red 
                rgb_cluster_array(index_array(i,2),index_array(i,1),1,index_array(i,3)) = 1; 
            elseif find(cluster_2_indices == i) 
                rgb_cluster_array(index_array(i,2),index_array(i,1),2,index_array(i,3)) = 1; 
            elseif find(cluster_3_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),3,index_array(i,3)) = 1;
            elseif find(cluster_4_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),1,index_array(i,3)) = 1;
                rgb_cluster_array(index_array(i,2),index_array(i,1),3,index_array(i,3)) = 1;
            elseif find(cluster_5_indices == i)
                rgb_cluster_array(index_array(i,2),index_array(i,1),2,index_array(i,3)) = 1;
                rgb_cluster_array(index_array(i,2),index_array(i,1),3,index_array(i,3)) = 1;
            end
        end
     rgb_stack = immovie(rgb_cluster_array); 
     implay(rgb_stack); 
     hold off; 
    
    end
    
%% Now we plot the clusters. 
figure; 

plot(clustering_variables(clusters==1,1),clustering_variables(clusters==1,2), ...
    'r.', 'MarkerSize', 12); 
hold on 
plot(clustering_variables(clusters==2,1),clustering_variables(clusters==2,2), ...
    'g.', 'MarkerSize', 12); 
hold on; 
plot(clustering_variables(clusters==3,1),clustering_variables(clusters==3,2), ...
    'b.', 'MarkerSize', 12); 
hold on; 
plot(clustering_variables(clusters==4,1),clustering_variables(clusters==4,2), ...
    'm.', 'MarkerSize', 12); 
hold on; 
plot(clustering_variables(clusters==5,1),clustering_variables(clusters==5,2), ...
    'c.', 'MarkerSize', 12); 

legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5'); 
xlabel('Peak Concenration [mM]'); 
ylabel('First Moment Integral, A.U.'); 
hold off; 

[fm_values,peak_values]=ginput(); 

min_fm = min(fm_values); 
max_fm = max(fm_values); 

min_peak = min(peak_values); 
max_peak = max(peak_values); 
% 

fm_indices = find(clustering_variables(:,1)<max_fm & clustering_variables(:,1)...
    >min_fm);

%[fm_indices] = find(clustering_variables(:,1),'max_peak','<','min_peak','>')

[peak_indices]=find(clustering_variables(:,2)<max_peak & clustering_variables(:,2)...
    >min_peak);

good_indices = intersect(fm_indices,peak_indices);

AIFs_2_average = zeros(dimt,length(good_indices)); 

%AIF_plot = zeros(dimx,dimy); 
for i = 1 : length(good_indices)
    % Find indices: 
    xind = index_array(good_indices(i),1); 
    yind = index_array(good_indices(i),2); 
    %zind = index_array(good_indices(i),3); 
    % Use indices to find appropriate CTCs
    if image_dims == 3 
        AIFs_2_average(:,i) = concentration_array(xind,yind,:); 
        rgb_cluster_array(yind,xind,:) = 1; 
    else 
        %find z indices
        zind = index_array(good_indices(i),3);
        AIFs_2_average(:,i) = concentration_array(xind,yind,zind,:); 
        rgb_cluster_array(yind,xind,:,zind) = 1; 
    end
      
end

%AIF_signal_2_average = zeros(length(time_vector), length(good_indices));
for i = 1:length(good_indices)
    %Find the signal values of the whole scan to plot
    xind = index_array(good_indices(i),1); 
    yind = index_array(good_indices(i),2); 
    if image_dims == 3
        AIF_signal_2_average(:,i) = image_array(xind, yind,:);
    else
        AIF_signal_2_average(:,i) = image_array(xind, yind, zind,:);
    end
end


 
%Now average the AIFs and plot a 'mean' AIF. 

mean_AIF = mean(AIFs_2_average,2); 
meanSignal = mean(AIF_signal_2_average,2);

if image_dims == 3
    %PLot the AIF
    figure; 
    plot(tx,mean_AIF);
    title('Mean AIF Plot');
    xlabel('Time'); 
    ylabel('[mM]'); 
    hold off; 
    % Now show the selected pixels 
    figure; 
    image(rgb_cluster_array);
    title('Selected Values in the Image Show in White')
    axis image; 
    hold off; 
    axis off;
    
    figure;
    time_vect_sec2 = 0 : length(meanSignal) -1; 
    time_vect_sec2 = time_vect_sec2 * TR;
    plot(time_vect_sec2 .* 60, meanSignal);
    title('AIF Mean Signal from Start to End of Scan')
    xlabel('Time (s)')
    ylabel('Signal Intensity (au)')
    
else
    %Plot the AIF
    figure; 
    plot(tx,mean_AIF);
    xlabel('Time'); 
    ylabel('[mM]'); 
    hold off; 
    figure;
    % Now display the 'selected voxels'
    figure; %new
     rgb_stack = immovie(rgb_cluster_array); 
     implay(rgb_stack); 
     hold off; 
     
     figure;
    plot(time_vector * 60, Signal_2_average);
    
     figure;
    time_vect_sec2 = 0 : length(meanSignal) -1; 
    time_vect_sec2 = time_vect_sec2 * TR;
    plot(time_vect_sec2 .* 60, meanSignal);
    title('AIF Mean Signal from Start to End of Scan')
    xlabel('Time (s)')
    ylabel('Signal Intensity (au)')
    
end

% Ask user if they are OK with it...

good_aif = menu('Are you OK with chosen AIF?', 'Yes', 'No');

end
%might need to transpose this final matrix 
meanAIF = mean_AIF;



    




 






