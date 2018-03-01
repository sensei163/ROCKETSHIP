function [CBF,CBV,MTT] = DSC_convolution_sSVD(concentration_array, AIF,deltaT,Kh,rho,Psvd,method,image_path) 

% The functioin DSC_CONVOLUTION is used to to calculate the residual function, R(t), and scaling factor, F,
%to be used in subsequent blood
% flow calculations.  This function uses a block circulant deconvolution
% method to solve for the residual function algebraically.  

% FUNCTION OUTPUTS: 

%   CBF: Cerebral blood flow in units of ml/100g/min. Computed for every voxel. 

%   CBV: Cerebral blood volume. 

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

%   Psvd: A threshold value that determines which elements to remove from
%   the SVD matrices. Values that are very low can skew computations in a
%   bad way.  

%Revision history: 
% 07/19/2015: M. Ellis edited function to produce MTT and time-to-peak maps. 
%%%%%%%%%%%%%%%%%%%%****FUNCTION BODY BEGIN HERE****%%%%%%%%%%%%%%%%%%%%%%

%% ARRANGE INPUTS ARRAYS AND INITIALIZE OUTPUT ARRAYS: 

% img = load_nii('~/desktop/DSC_test_image.nii'); 
% concentration_array = img.img; 
% concentration_array(concentration_array<0) = 0; 
% concentration_array(concentration_array >32767) = 0; 
% Kh = 0.71;   % difference in hematocrit between capillaries and large vessels. 
% rho = 1.04;  % Tissue density. In g/ml. But we want 100g/ml. So we divide by 100;  
% AIF = mean_AIF; 
% deltaT = TR/60;  % We want minutes.  
% Psvd = 0.001; 

if ndims(concentration_array) ==4
 [dimx, dimy, dimz, dimt] = size(concentration_array); 
else 
 [dimx,dimy,dimt] = size(concentration_array); 
 dimz = 1; 
end
 
 CBF = zeros(dimx,dimy,dimz); 
 CBV = zeros(dimx,dimy,dimz); 
 MTT = zeros(dimx,dimy,dimz); 
 TTP = zeros(dimx,dimy,dimz); 
 
 if dimz == 1 
     CBF = squeeze(CBF); 
     CBV = squeeze(CBV);
     MTT = squeeze(MTT); 
     TTP = squeeze(TTP); 
     
 end 
 

 
 % We need AIF to be in the form of a lower trianguar matrix, A, in order
 % to perform our subsequent singular value decomposition. Also, we prefilter A and account for discretization errors. 
 % This is performed below.
 
 A = zeros(dimt); 
 
 for i = 1 : dimt
     
     for j = 1 : dimt 
         
         if i>j && 1<i && i<dimt && 1<j && j<dimt 
             
         A(i,j) = deltaT*(AIF(i-j)+4*AIF(i-j+1)+AIF(i-j+2))/6;  
         
         elseif i == j 
             
         A(i,j) = deltaT*(AIF(i-j+1)); 
           
         end 
     end
 end
 
 % We are also concerned with the CBV. To calculate the CBV we need to
 % integrate the AIF: 
 
 % Defining a time vector to integrate according to: 
 time_vect = 0 : deltaT : (dimt-1)*deltaT; 
 time_vect = time_vect'; 
 AIF_int = trapz(time_vect,AIF);
 
%% SOLVE FOR THE RESIDUAL ON A VOXEL-BY-VOXEL BASIS: 

            % STANDARD SINGULAR VALUE DECOMPOSITION  
            % We are trying to solve a linear system that assume the form: 
            % c = A * b; 
            % Where: c = Ct; A = AIF matrix; b = Rt * Ft; 
            
            % We decompose A, take its inverse, and solve for b. 
            
            [U, S, V] = svd(A); %(2:23,2:23)); 
            S(S<Psvd) = 0; 
for k = 1 : dimz 
    for j = 1 : dimy 
        for i = 1 : dimx 
            
            % A concentration array: 
            if ndims(concentration_array) == 4
            c = concentration_array(i,j,k,:); 
            c = squeeze(c); 
            else 
            c = concentration_array(i,j,:);
            c = squeeze(c); 
            end
            % Now we decompose D using singular value decomposition: 
         
            % The above uses the following equality: D = Uc * Sc * Vc'; 
            
            % With the decomposed values we can solve for the residual
            % function, R(t); 
            
            % 1: Define W. W = 1/Sc along its diagonal and zeros otherwise. 
            % 2: Take the transpose of Uc  
            
            W = diag(diag(1./S));   % Step 1. 
            W(W==Inf) = 0;
            
            UT = U';                % Step 2. 
            
            % We now have the matrices needed to solve our system
            % algebrically: 
            
            f = V * W * (UT * c);        
            
            %Commented out 07/19/2015; This is solved for below. 
            %CBF(i,j,k) =  Kh/rho .* max(max(f));   
            
            % To calculate the CBV we need to integrate c; 
            c_int = trapz(time_vect,c); 
            
            % 07/19/2015: Added some code to produce time to peak maps
            if method ==1 
                if ndims(concentration_array) == 4
                  CBF(i,j,k) =  Kh/rho .* max(max(f));                      % 07/19/2015: Will this work for the 4-D case? I'm not sure...
                  CBV(i,j,k) = (Kh/rho) * (c_int / AIF_int); 
                  MTT(i,j,k) = CBV(i,j,k) / CBF(i,j,k); 
                  [~,max_ind] = max(squeeze(concentration_array(i,j,k,:))); 
                  TTP(i,j,k) = time_vect(max_ind); 
                  %We need to threshold the TTP to eliminate noise voxels. Will work with Axel on this one soon... 
                else 
                  CBF(i,j) =  Kh/rho .* max(max(f));
                  CBV(i,j) = (Kh/rho) * (c_int / AIF_int); 
                  MTT(i,j) = CBV(i,j) / CBF(i,j); 
                  [~,max_ind] = max(squeeze(concentration_array(i,j,:))); 
                  TTP(i,j) = time_vect(max_ind);
                end
            elseif method == 2 
                if ndims(concentration_array) == 4
                  CBF(i,j,k) =  Kh/rho .* max(max(f)); 
                  CBV(i,j,k) = (Kh/rho) .* max(max(f)) .* trapz(time_vect,f); 
                  MTT(i,j,k) = trapz(time_vect,f); 
                  [~,max_ind] = max(squeeze(concentration_array(i,j,k,:)));
                  TTP(i,j,k) = time_vect(max_ind); 
                else 
                  CBF(i,j) =  Kh/rho .* max(max(f));
                  CBV(i,j) = (Kh/rho) .* max(max(f)) .* trapz(time_vect,f); 
                  MTT(i,j) = trapz(time_vect,f); 
                  [~,max_ind] = max(squeeze(concentration_array(i,j,:))); 
                  TTP(i,j) = time_vect(max_ind);
                end
            end  
        end
    end
end

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

    
     


    

