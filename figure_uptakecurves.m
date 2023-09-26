
clear all; clf
addpath Uptake
fileList = {'sm006779.mat','sm006784.mat','sm006785.mat','sm006786.mat'};
temp_list = {'120 K', '130 K', '135 K', '140 K','150 K','160 K'};
pre_t_list = {480, 230, 280,150, 45, 68}; %A time before water deposition starts

startExp = {421,204,254,124}; %Here water deposition starts
figure=figure('units','centimeters','position',[4,2,16,14],'color','white','DefaultLineLineWidth',2,'DefaultTextFontSize',18);

for i=1:length(fileList)
    load(fileList{i});
    pre_t = pre_t_list{i};
    
    startIndex = cell2mat(startExp(i)); 
    
    pressure = meas.samp_pressure_e; %pressure extraction gauge  
    pre_press = mean(pressure(meas.t<pre_t)); %average pressure before deposition 

    temperature = mean(meas.samp_temp(meas.t>pre_t));

    dexposure = ((pressure(1:end-1)+pressure(2:end))/2-pre_press).*(meas.t(2:end)-meas.t(1:end-1)); 
    
    exposure = cumsum(dexposure);
    exposure_langmuir=exposure*10^6;

    signal = (meas.counts(1:end-1)+meas.counts(2:end))/2;
    
    %Normalise and plot from startExp
    signal = signal./signal(startIndex);
    semilogy(exposure_langmuir(startIndex:end),signal(startIndex:end),'displayname', temp_list{i})
    xlim([0 0.12])
    ylabel('Normalized counts')
    xlabel('H_2O exposure (L)')


    hold on
    
   
end
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
legend
grid
hold off


