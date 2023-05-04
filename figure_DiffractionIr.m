
clear all; clf
fileList = {'as010376.mat','as010377.mat'}
names = {'Ir(111)','H_2O Ir(111)'};
colourList = 'rbgcmk';
figure=figure('units','centimeters','position',[4,2,16,14],'color','white','DefaultLineLineWidth',2,'DefaultTextFontSize',18);

addpath Diffraction

for i=1:length(fileList)
    load(fileList{i})
    signal = (meas.scan.counts(1:end-1)+meas.scan.counts(2:end))/2;
    signal = signal./max(signal);
    name = names{i};
    semilogy(meas.scan.angles-meas.ga,meas.scan.counts,colourList(i),'displayname', [num2str(name)]) %Plot the angles as relative to specular 

    hold on 
end
xlim([22-meas.ga 53-meas.ga])

xlabel('\Delta\theta_i relative to specular (\circ)')
ylabel('Scattered signal (arb units)')
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
legend
grid
hold off


