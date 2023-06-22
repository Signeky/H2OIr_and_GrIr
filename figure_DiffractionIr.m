
clear all; close all; 
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
    semilogy(meas.scan.angles,meas.scan.counts,colourList(i),'displayname', [num2str(name)]) %Plot the angles as relative to specular 

    hold on 
end
xlim([22 53])

xlabel('Angle of incidence, \theta_i (\circ)')
ylabel('Scattered signal (arb units)')
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
legend
grid
hold off



%% Plot with Delta K rather than angle

% Calibration of the beam 08/12/2021 12:19
load('beam_20211208_1219.mat')
beam.E0
% Calibration of the beam 15/12/2021 16:08
load('beam_20211215_1608.mat')
beam.E0
% The 2 diffraction scans being plotted in this script were taken 
% 9/12/2021, 13:08 & 13:45
% Use the calibration from 08/12/2021
load('beam_20211208_1219.mat')
E0 = beam.E0; % meV

E0 = E0*1.603e-19/1e3; % J
h = 6.626e-34; % Planck constant, Js
m = 3.0160293*1.66054e-27; % Helium 3 rest mass, kg
k = (2*pi*sqrt(2*E0*m)/h); % m^-1, incident wavector
k = k/1e10; % A^-1
disp(['k = ' num2str(k) ' A^-1'])

load(fileList{1})
Ir_gamma = meas.scan.angles-meas.ga;
Ir_counts = meas.scan.counts;
Ir_DK = -2*k*cosd(22.2)*sind(Ir_gamma);

load(fileList{2})
IrH20_gamma = meas.scan.angles-meas.ga;
IrH20_counts = meas.scan.counts;
IrH20_DK = -2*k*cosd(22.2)*sind(IrH20_gamma);




semilogy(-Ir_DK, Ir_counts, colourList(1), 'displayname', names{1})
hold on
semilogy(-IrH20_DK, IrH20_counts, colourList(2), 'displayname', names{2})
hold off
xlabel('\Delta K (Å^{-1})')
ylabel('Scattered signal (arb units)')
xlim([min(IrH20_DK), 0])
legend
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
grid on
xlim([0.11 3.4])
ylim([2e-11 1e-8])

% Find the diffraction peak to measure the lattice constant
ind = Ir_DK < -2.5 & Ir_DK > -2.9;
y = flip(Ir_counts(ind));
x = flip(Ir_DK(ind));
pks = findpeaks(y, x);
G = abs(x(y == pks))
a = (2*pi/G)/cosd(30) % Accound for measuring along the Gamma-M direction
disp(['Measured lattice constant of Ir(111) = ' num2str(a) ' Å'])
% Note spacing of data points is ~0.02 A^(-1), giving ~0.02 A
% uncertainty in the lattice constant. The beam energy also drifts slightly
% 0.02meV over a week, so this is a significantly smaller uncertainty.
% I think we can quote to 3 s.f. FWHM of the diffraction peak is ~0.07A^-1
% which would give a standard deviation of the peak (assuming Gaussian) of
% ~0.03A^-1. Translating to ~0.03A in the lattice constant.
