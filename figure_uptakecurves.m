% use the following to determine the pre-time
%
% clf; load('sm006793.mat'); semilogy(meas.t,meas.samp_pressure_e,meas.t,meas.counts) 

clear all; clf
addpath Uptake
%pressure measured from extraction gauge 145 K file is corrupt and does not show.  
fileList = {'sm006779.mat','sm006784.mat','sm006785.mat','sm006786.mat'};
temp_list = {'120 K', '130 K', '135 K', '140 K','150 K','160 K'};
pre_t_list = {480, 230, 280,150, 45, 68};
%max_p_list = {'1.1E-8 cst','1.1E-8 cst','1.2E-8 spike?','1.4E-8 cst ','3.3E-8 slw','2E-8 cst','','112 K','','','',''};

startExp = {421,204,254,124};
figure=figure('units','centimeters','position',[4,2,16,14],'color','white','DefaultLineLineWidth',2,'DefaultTextFontSize',18);

for i=1:length(fileList)
    load(fileList{i});
    pre_t = pre_t_list{i};
    
    startIndex = cell2mat(startExp(i)); 
    
    pressure = meas.samp_pressure_e; %pressure extraction gauge    
    pre_press = mean(pressure(meas.t<pre_t)); %average pressure before deposition 

    temperature = mean(meas.samp_temp(meas.t>pre_t));
    
    dexposure = ((pressure(1:end-1)+pressure(2:end))/2-pre_press).*(meas.t(2:end)-meas.t(1:end-1)); %-pre_press causes disortoion
    
    
    exposure = cumsum(dexposure);
    
    signal = (meas.counts(1:end-1)+meas.counts(2:end))/2;
    
    %Normalise and plot from startExp
    signal = signal./signal(startIndex);
    semilogy(exposure(startIndex:end),signal(startIndex:end),'displayname', temp_list{i})
    xlim([0 1.2*10^(-7)])
    ylabel('Normalized counts (arb units)')
    xlabel('H_2O exposure (mbar x s)')


    hold on
    
   
end
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
legend
grid
hold off


