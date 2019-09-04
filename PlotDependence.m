%% Plot SPR and distance dependent/power depedent absolute photon count rate. 

directories = dir();

dirFlags = [directories.isdir];
subFolders = directories(dirFlags);

pcles = length(subFolders) - 2; 

for i = 3 : length(subFolders)
    cd(subFolders(i).name)
    load PCR
    load decayrates
    NonenhDistDepend_HI{i - 2} = PCR.NonenhCountRateDistQYdepend_HI;
    EnhDistDepend_HI{i -2} = PCR.EnhCountRateDistQYdepend_HI;
    NonenhDistDepend_LO{i - 2} = PCR.NonenhCountRateDistQYdepend_LO;
    EnhDistDepend_LO{i -2} = PCR.EnhCountRateDistQYdepend_LO;
    NonenhPowerQYDepend{i - 2} = PCR.NonenhCountRatePowerQYdepend;
    EnhPowerQYDepend{i - 2} = PCR.EnhCountRatePowerQYdepend;
    I_satenh(i - 2) = PCR.I_satenh(11)./decayrates.ee(11);
    maxPCRenh (i - 2) = PCR.maxPCRenh(6);
    maxPCRenhfactor(i - 2) = PCR.maximumPCRenhancement(6); 
    WL(i - 2) = 1248/decayrates.Lorentz(3);
    cd ..
end 

%% Plot SPR dependent, power dependent, QY = 0.05
figure
x = PCR.I_inc ;
pcles_id = [ 3, 6, 9, 11];
for i = 1 : length(pcles_id)
    loglog(x, EnhPowerQYDepend{pcles_id(i)}(6,:))
    %vline(I_satenh(pcles_id(i)))
    %hline(maxPCRenh(pcles_id(i)))
    hold on 
end
nonenh_PCR = PCR.NonenhCountRatePowerQYdepend(6,:) ;
loglog(x, nonenh_PCR,'-')
%vline(PCR.I_sat)
%hline(PCR.maxPCR(6))
xlabel('Power Wm^{-2}')
ylabel('Photon count rate')
lg = split(num2str(round(WL(pcles_id))),'  ');
for i = 1 : length(lg)
    lg{i} = [lg{i},' nm'];
end 
lg{end + 1} = 'Non-enhanced';
legend( lg )
MinWhitSpace
saveas( gcf,'AbsPCR_power_depend.fig' )
saveas( gcf,'AbsPCR_power_depend.png' )
save incidentpowerdensity x 
save nonenh_PCR_QY0.5.mat nonenh_PCR
%% Plot SPR dependent, at fixed powers, QY = 0.05
figure
idx = [36, 41,43, 46, 48, 50];
power = PCR.I_inc(idx);
for i = 1 : pcles
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

%% Plot SPR dependent brightness enhancement, QY = 0.05
figure
plot(WL,maxPCRenhfactor,'o-')
xlabel( 'SPR wavelength (nm)' )
ylabel('Brightness enhancement')
save WL WL
save maxPCRenhfactor maxPCRenhfactor
saveas( gcf,'SPR_vs_Brightness_enh.fig' )
saveas( gcf,'SPR_vs_Brightness_enh.png' )


%% Plot SPR dependent, at fixed powers, QY = 1

figure
idx = [36, 41,43, 46, 48, 50];
power = PCR.I_inc(idx);
for i = 1 : pcles
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

% plot SPR of maximum enhancement as a function power
figure
[maxi, id_max_1 ] = max( enhpqydepend(:,1) );
[maxi, id_max_2 ] = max( enhpqydepend(:,2) );
[maxi, id_max_3 ] = max( enhpqydepend(:,3) );
[maxi, id_max_4 ] = max( enhpqydepend(:,4) );
[maxi, id_max_5 ] = max( enhpqydepend(:,5) );
[maxi, id_max_6 ] = max( enhpqydepend(:,6) );
id_max = [WL(id_max_1),WL(id_max_2),WL(id_max_3),WL(id_max_4),WL(id_max_5),WL(id_max_6)];
plot(PCR.I_inc(idx), id_max )

save id_max id_max
save idx idx

%% Plot SPR and distance dependent SPR. Distance = 1, 3, 10, 30, 100. 

idx_d = [ 3, 10, 35, 90, 200 ];

figure
for i = 1 : pcles
    loglog(decayrates.d_BEM, EnhDistDepend_HI{i}(7,:))
    hold on 
end 
legend(lg)
hline(NonenhDistDepend_HI{1,7}(7),'r:','Non-enhanced' )
xlabel( 'd (nm)' )
ylabel('Photon count rate')

MinWhitSpace
saveas( gcf,'Distance_depend.fig' )
saveas( gcf,'Distance_depend.png' )

%% Plot 2D correlation of SPR WLs (x) - Distance (y) - PCR (z) heat maps, fixed QY = 0.5

% At high power

figure
for i = 1 : pcles
    z(i,:) = EnhDistDepend_HI{i}(6,:);
end 
    
    
imagesc( WL, decayrates.d_BEM, z');
colorbar
colormap jet
set(gca,'ColorScale','log')
xlabel('SPR wavelength (nm)')
ylabel('molecule-particle separation d (nm)')
zlabel('Photon s^{-1}')

saveas(gcf,'SPR-d-PCR correlation HIGH exc.fig')
saveas(gcf,'SPR-d-PCR correlation HIGH exc.png')

% At low power

figure
for i = 1 : pcles
    z(i,:) = EnhDistDepend_LO{i}(6,:);
end 
    
    
imagesc( WL, decayrates.d_BEM, z');
colorbar
colormap jet
set(gca,'ColorScale','log')
xlabel('SPR wavelength (nm)')
ylabel('molecule-particle separation d (nm)')
zlabel('Photon s^{-1}')

saveas(gcf,'SPR-d-PCR correlation LOW exc.fig')
saveas(gcf,'SPR-d-PCR correlation LOW exc.png')

%% Plot distance-dependent brightness enhancement, at HIGH power
figure
for i = 1 : pcles
    z(i,:) = EnhDistDepend_HI{i}(6,:)./NonenhPowerQYDepend{1,1}(6,end);
end 

imagesc( WL, decayrates.d_BEM, z');
colorbar
colormap jet
set(gca,'ColorScale','log')
xlabel('SPR wavelength (nm)')
ylabel('molecule-particle separation d (nm)')
zlabel('Brightness enhancement')

saveas(gcf,'SPR-d-PCR correlation enhancement HIGH exc.fig')
saveas(gcf,'SPR-d-PCR correlation enhancement HIGH exc.png')