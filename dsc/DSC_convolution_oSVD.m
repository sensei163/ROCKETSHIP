function [CBF,CBV,MTT] = DSC_convolution_oSVD(concentration_array, AIF,deltaT,Kh,rho,OI,method,image_path) 

% The functioin DSC_CONVOLUTION is used to to calculate the residual function, R(t), and scaling factor, F,
%to be used in subsequent blood
% flow calculations.  This function uses a block circulant deconvolution
% method to solve for the residual function algebraically.  

% FUNCTION OUTPUTS: 

%   CBF: Cerebral blood flow in units of ml/(100g * min). Computed for every voxel. 

%   CBV: Cerebral blood volume in units of ml/(100g * min). 

%   MTT: Mean transit time of tracer. 

% FUNCTION INPUTS: 

%   CONCENTRATION_ARRAY: The array of concentration values. These are C(t)
%   values that are used to compute the residual. 

%   AIF: Is the arterial input function. It is concentration values of
%   arterials. This tells use how much contrast agent in flowing into the
%   tissue region.  

%   DeltaT: The time steps between scan points. We are essentially performing
%   a numerical integration to calculate cerebral blood flow and cerebral
%   blood volume. 

%   Kh: A constant ratio of capillary to artery hematocrit

%   RHO: Average tissue density.  (Its very close to unity and should
%   remain there unless you know what you're doing. 

%   OI: A threshold value that determines the degree of oscilliations allowed
%   in the final residual function. The Psvd is increased until the OI 

%   method: 1 calculates CBV as CBV= 100*(Kh/rho) * (c_int / AIF_int)
%   2 calculates CBV as CBV= 100*(Kh/rho) .* (integral of the residual fxn over the course of the scan)

%Revision history: 
% 07/19/2015: M. Ellis edited function to produce MTT and time-to-peak maps. 
%%%%%%%%%%%%%%%%%%%%****FUNCTION BODY BEGIN HERE****%%%%%%%%%%%%%%%%%%%%%%

%% ARRANGE INPUTS ARRAYS AND INITIALIZE OUTPUT ARRAYS:
if ndims(concentration_array) == 4
 [dimx, dimy, dimz, dimt] = size(concentration_array); 
else 
 [dimx,dimy,dimt] = size(concentration_array); 
 dimz = 1; 
end
 
 CBF = zeros(dimx,dimy,dimz); 
 CBF_int = zeros(dimx,dimy,dimz);
 CBV = zeros(dimx,dimy,dimz); 
 MTT = zeros(dimx,dimy,dimz); 
 TTP = zeros(dimx,dimy,dimz); 
 
 if dimz == 1 
     CBF = squeeze(CBF); 
     CBF_int = zeros(dimx,dimy,dimz);
     CBV = squeeze(CBV);
     MTT = squeeze(MTT); 
     TTP = squeeze(TTP); 
     
 end 
 
 % We are also concerned with the CBV. To calculate the CBV we need to
 % integrate the AIF: 
 
 % Defining a time vector to integrate according to: 
 time_vect = 0 : deltaT : (dimt-1)*deltaT; 
 time_vect = time_vect'; 
 
 %Define a time vector that extends the length of f (double the amount of
 %time as time_vect plus the time resolution of the MRI scan)
 time_vect_ext = 0 : time_vect(end)/(numel(time_vect)-1) : (time_vect(end) * 2) + time_vect(2);
 time_vect_ext = time_vect_ext';
 AIF_int = trapz(time_vect,AIF);

 %% Linearize AIF time curve
 %linarization function rewriten 5/23/18 Prince to match Ostergaard et al.,1996

 
 %zero pad AIF 
 AIF(2*dimt) = 0; 
 
 A = zeros(2*dimt);
 
 for i = 1 : 2*dimt
     for j = 1 : 2*dimt
         if i == j 
             A(i,j) = (deltaT*(4*AIF(i-j+1)+AIF(i-j+2)))/6; 
             %accounts for when the AIF(i-j) term equals AIF(0) = 0 
               %AIF(0) refers to the time before the injection when the
               %concentration should be zero 
         elseif j < i && i == 2*dimt && j == 1
             A(i,j) = (deltaT*(AIF(i-j)+4*AIF(i-j+1)+AIF(2*dimt)))/6;
             %i-j exceeds indexing assume AIF(dimt+1) = AIF(dimt) 
         elseif j < i
              A(i,j) = (deltaT*(AIF(i-j)+4*AIF(i-j+1)+AIF(i-j+2)))/6;
         end 
     end
 end
 
 L = 2*dimt;
 D = zeros(L);
 
 for i = 1 : L
     for  j = 1 : L
         if j <= i
             D(i,j) = A(i,j);
         else
             D(i,j) = A(L + i - j + 1,1); 
         end
     end    
 end
 
 
 
%% SOLVE FOR THE RESIDUAL ON A VOXEL-BY-VOXEL BASIS AND CHECK FOR AGREEMENT WITH THE OI (rerun with AN UPDATED Psvd IF NECESSARY):
            
for k = 1 : dimz 
    for j = 1 : dimy 
        for i = 1 : dimx 
            
            %create a psuedo-do/while loop to increment Psvd
            threshold_met = false; Psvd_temp = 0; %what percentage cut-off the program starts with 
           
            %oSVD finds the lowest Psvd that gives produces a residual function with an osscilatory index
            %below the desired value (OI). To do this Psvd is increased
            %incrementally starting from this value. 

            while(threshold_met == false)
            
            % STANDARD SINGULAR VALUE DECOMPOSITION  
            % We are trying to solve a linear system that assumes the form: 
            % c = A * b; 
            % Where: c = Ct; A = AIF matrix; b = Rt * Ft; 
            
            % We decompose A, take its inverse, and solve for b. 
            
            [U, S, V] = svd(D); %(2:23,2:23)); 
            Psvd_threshold = max(max(max(S))) * (Psvd_temp/100); 
            S(S<Psvd_threshold) = 0; %thresholding of diaganal values to eliminate the effect of noise
            
            % A concentration array: 
            if ndims(concentration_array) == 4
            c = concentration_array(i,j,k,:); 
            c = squeeze(c); 
            else 
            c = concentration_array(i,j,:);
            c = squeeze(c); 
            end
            
            %zero pad c 
            c(2*dimt) = 0; 
            
            % Now we decompose D using singular value decomposition: 
         
            % The above uses the following equality: D = Uc * Sc * Vc'; 
            
            % With the decomposed values we can solve for the residual
            % function, R(t); 
            
            % 1: Define W. W = 1/Sc along its diagonal and zeros otherwise. 
            % 2: Take the transpose of Uc  
            
            W = diag(diag(1./S));   % Step 1. 
            W(isinf(W)) = 0;
            
            UT = U';                % Step 2. 
            
            % We now have the matrices needed to solve our system
            % algebrically: 
            
            f = V * W * (UT * c);        
           
            % To calculate the CBV we need to integrate c; 
            c_int = trapz(time_vect,c(1:dimt)); 
            
            %test r for threshold, rerun with different Psvd or calculate
            %maps (from Wu et al., 2003)
            
            sum = 0;
            for sum_iterator = 2:(L-1)
                sum = sum + abs(f(sum_iterator +1) - 2*f(sum_iterator) + f(sum_iterator-1));
            end
            O = (1/L)*(1/max(max(f)))*sum; 
            
            %screen for threshold and for NaN/Inf, the latter case
            %indicates a max f of zero and most likely a background pixel,
            %for which the O value doesn't matter
            if O > OI && (~isnan(O)) && (~isinf(O))
                Psvd_temp = Psvd_temp + 2.5;
            else
            
                threshold_met = true;
               
                res_array(:,i,j) = f;
                
                %to calcualte CBF we need the max of R(t) and the integral of R(t)
                %from zero to infinity 
                max_r = max(f); %handles 3D and 4D case max(b) should be suffienct
                r_int = trapz(time_vect,f(1:numel(time_vect))); 
                
            % 07/19/2015: Added some code to produce time to peak maps
          if method ==1 
                if ndims(concentration_array) == 4
                  CBF(i,j,k) =  100*(Kh/rho)*max(max(f)); %Knuttson et al. and Fieselmann et al., 2011
                  CBV(i,j,k) = 100*(Kh/rho) * (c_int / AIF_int); %Knuttson et al. 
                  MTT(i,j,k) = (CBV(i,j,k)) / CBF(i,j,k); %Knutsson et al. (*100 will cancel)
                  [~,max_ind] = max(squeeze(concentration_array(i,j,k,:))); %Fieselmann et al., 2011
                  TTP(i,j,k) = time_vect(max_ind); %in mins %Fieselmann et al., 2011
                  %We need to threshold the TTP to eliminate noise voxels. Will work with Axel on this one soon... 
                else 
                  CBF(i,j) =  100*(Kh/rho)*max(max(f)); %Knuttson et al. and Fieselmann et al., 2011
                  CBV(i,j) = 100*(Kh/rho) * (c_int / AIF_int); %Knutsson et al. 
                  CBV_test(i,j) = 100*(Kh/rho) .* trapz(time_vect_ext,f);
                  MTT(i,j) = (CBV(i,j)) / CBF(i,j); %Knutsson et al. (*100 in both will cancel)
                  [~,max_ind] = max(squeeze(concentration_array(i,j,:))); %Fieselmann et al., 2011
                  TTP(i,j) = time_vect(max_ind); %in mins %Fieselmann et al., 2011
                end
            elseif method == 2 
                if ndims(concentration_array) == 4
                  CBF(i,j,k) = 100*(Kh/rho)* max(max(f)); %Knutsson et al., 2010 and %Fieselmann et al., 2011
                  %CBV(i,j,k) = (Kh/rho) .* max(max(f)) .* trapz(time_vect,f); 
                  CBV(i,j,k) = 100*(Kh/rho) .* trapz(time_vect_ext,f); %Fieselmann et al., 2011
                  MTT(i,j,k) = trapz(time_vect,f)/max(max(f)); %Fieselmann et al., 2011
                  [~,max_ind] = max(squeeze(concentration_array(i,j,k,:))); %Fieselmann et al., 2011
                  TTP(i,j,k) = time_vect(max_ind); %in mins %Fieselmann et al., 2011
                else 
                  CBF(i,j) = 100*(Kh/rho)*max(max(f)); % Knutsson et al., 2010 and Fieselmann et al., 2011
                  %CBV(i,j) = (Kh/rho) .* max(max(f)) .* trapz(time_vect,f); 
                  CBV(i,j) = 100*(Kh/rho) .* trapz(time_vect_ext,f); %Fieselmann et al., 2011
                  MTT(i,j) = trapz(time_vect,f)/max(max(f)); %Fieselmann et al., 2011 
                  [~,max_ind] = max(squeeze(concentration_array(i,j,:))); %Fieselmann et al., 2011
                  TTP(i,j) = time_vect(max_ind); %in mins %Fieselmann et al., 2011
                end  
            end
        end
    end
end

end

%save the generated maps

CBF_map = make_nii(CBF); 
cbf_file = strcat(image_path,'CBF_map.nii'); 
save_nii(CBF_map, cbf_file); 

CBV_map = make_nii(CBV); 
cbv_file = strcat(image_path,'CBV_map.nii'); 
save_nii(CBV_map, cbv_file); 
 
MTT_map = make_nii(MTT); 
mtt_file = strcat(image_path,'MTT_map.nii'); 
save_nii(MTT_map, mtt_file); 

TTP_map = make_nii(TTP); 
ttp_file = strcat(image_path,'TTP_map.nii'); 
save_nii(TTP_map, ttp_file); 

    
end


    



