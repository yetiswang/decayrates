%% Plot SPR and distance dependent/power depedent absolute photon count rate. 

directories = dir();

dirFlags = [directories.isdir];
subFolders = directories(dirFlags);

for i = 3 : 14
    cd(subFolders(i).name)
    load PCR
    NonenhDistDepend_HI{i - 2} = PCR.NonenhCountRateDistQYdepend_HI;
    EnhDistDepend_HI{i -2} = PCR.EnhCountRateDistQYdepend_HI;
    NonenhDistDepend_LO{i - 2} = PCR.NonenhCountRateDistQYdepend_LO;
    EnhDistDepend_LO{i -2} = PCR.EnhCountRateDistQYdepend_LO;
    NonenhPowerQYDepend{i - 2} = PCR.NonenhCountRatePowerQYdepend;
    EnhPowerQYDepend{i - 2} = PCR.EnhCountRatePowerQYdepend;
    cd ..
end 

%% Plot SPR dependent, power dependent, QY = 0.5
figure
x = PCR.I_inc ;
for i = 1 : 12
    plot(x, EnhPowerQYDepend{i}(6,:))
    hold on 
end 
plot(x, PCR.NonenhCountRatePowerQYdepend(6,:),':')
xlabel('Power Wm^{-2}')
ylabel('Photon count rate')
load WL
lg = split(num2str(round(WL)),'  ');
lg{end + 1} = 'Non-enhanced';
legend( lg )
MinWhitSpace
saveas( gcf,'AbsPCR_power_depend.fig' )
saveas( gcf,'AbsPCR_power_depend.png' )
%% Plot SPR dependent, at fixed powers, QY = 0.5
figure
idx = [36, 41,43, 46, 48, 50];
power = PCR.I_inc(idx);
for i = 1 : 12
    for j = 1 : length(power)
        enhpqydepend(i,j) = EnhPowerQYDepend{i}(6,idx(j));
    end
end
str = cell(1, length(idx));
for j = 1 : length(idx)
    semilogy(WL,enhpqydepend(:,j),'o-')
    hold on
    str{j} = num2str(PCR.I_inc(idx(j)),'%.E');
end 
legend(strcat(str, ' W m^{-2}'))
xlabel( 'SPR wavelength (nm)' )
ylabel('Photon count rate')

saveas( gcf,'AbsPCR_SPR_depend.fig' )
saveas( gcf,'AbsPCR_SPR_depend.png' )

%% Plot SPR dependent, at fixed powers, QY = 1

figure
idx = [36, 41,43, 46, 48, 50];
power = PCR.I_inc(idx);
for i = 1 : 12
    for j = 1 : length(power)
        enhpqydepend(i,j) = EnhPowerQYDepend{i}(7,idx(j));
    end
end
str = cell(1, length(idx));
for j = 1 : length(idx)
    semilogy(WL,enhpqydepend(:,j),'o-')
    hold on
    str{j} = num2str(PCR.I_inc(idx(j)),'%.E');
end 
legend(strcat(str, ' W m^{-2}'))
xlabel( 'SPR wavelength (nm)' )
ylabel('Photon count rate')

saveas( gcf,'AbsPCR_SPR_depend.fig' )
saveas( gcf,'AbsPCR_SPR_depend.png' )

%% Plot SPR and distance dependent SPR. Distance = 1, 3, 10, 30, 100. 

idx_d = [ 3, 10, 35, 90, 200 ];

figure
for i = 1 : 12
    loglog(d_BEM, EnhDistDepend_HI{i}(6,:))
    hold on 
end 
legend(lg)
hline(NonenhDistDepend_HI{1,6}(6),'r:','Non-enhanced' )
xlabel( 'd (nm)' )
ylabel('Photon count rate')

MinWhitSpace
saveas( gcf,'Distance_depend.fig' )
saveas( gcf,'Distance_depend.png' )

%% Plot 2D correlation of SPR WLs (x) - Distance (y) - PCR (z) heat maps, fixed QY = 0.5

% At high power

figure
for i = 1 : 12
    z(i,:) = EnhDistDepend_HI{i}(6,:);
end 
    
    
imagesc( WL, d_BEM, z');
colorbar
colormap jet
xlabel('SPR wavelength (nm)')
ylabel('molecule-particle separation d (nm)')
zlabel('Photon s^{-1}')

saveas(gcf,'SPR-d-PCR correlation HIGH exc.fig')
saveas(gcf,'SPR-d-PCR correlation HIGH exc.png')

% At low power

figure
for i = 1 : 12
    z(i,:) = EnhDistDepend_LO{i}(6,:);
end 
    
    
imagesc( WL, d_BEM, z');
colorbar
colormap jet
xlabel('SPR wavelength (nm)')
ylabel('molecule-particle separation d (nm)')
zlabel('Photon s^{-1}')

saveas(gcf,'SPR-d-PCR correlation LOW exc.fig')
saveas(gcf,'SPR-d-PCR correlation LOW exc.png')

