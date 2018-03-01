# DSC_GUI

AIF_input function formating: 

AIF Input File Formate:
  Contains 1 Arrays: meanAIF
    
    meanAIF: an array of the average AIF values
    
    bolus_time the a slice number or index with the bolus is injected

    bolus_time = contrast injection index - start index of AIF 
    
    ex: if the contrast is injected at the same time as the AIF then the 
    bolus_time is zero 

        if AIF starts at time zero then bolus_time = start index of the AIF
        
        if the AIF starts at a nonzero index then use the formula 
        bolus_time = contrast injection index - start index of AIF to get
        the appropriate bolus_time value 
