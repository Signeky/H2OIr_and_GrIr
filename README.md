## HeSE data for water on iridium(111) and Gr/iridium(111)
This project includes raw helium-3 spin-echo data and the analysis of the adsorption and dynamics of water on iridium(111) and graphene/iridium(111). MatLab R2021a Update 4 is used for the analysis. 


###Dynamics measurements

Plot_dynamics_Ir125.m MATLAB 	script to load and analyse the raw polarised data collected on Ir(111) (GM azimuth, 125 K)
Plot_dynamics_Ir135.m MATLAB 	script to load and analyse the raw polarised data on Ir(111)  (GM azimuth, 135 K)
Plot_dynamics_GrIr125.m MATLAB 	script to load and analyse the raw polarised data collected on GrIr(111) (GM azimuth, 125 K)
  
figure_dynamics.m MATLAB script to load and analyse the raw polarisation data for all substrates and produce a dynamic figure 
figure_Arrhenius.m MATLAB script to load raw polarisation on Ir(111) data and calculate the Arrhenius values
	
/Dynamics/Ir_and_GrIr	folder containing polarisation curves of the raw spin-echo measurements for the Ir(111) and GrIr(111) substrates
/Dynamics/GrNi		data from  DOI https://doi.org/10.1038/s41467-021-23226-5

    
    
Water uptake curves
	
	figure_uptakecurves.m MATLAB script to load and analyse the raw data of the measured uptake curves and produce the uptake figure
	
	/Uptake/ raw data of uptake curves 


diffraction measurements 

	figure_DiffractionIr.m MATLAB script to load and analyse the diffraction scan on Ir(111) 
	figure_DiffractionGrIr.m MATLAB script to load and analyse the diffraction scan on GrIr(111) 
	
	/Diffraction/ raw data of diffraction scans


Information on the measurement conditions is included in the MATLAB data files. 
