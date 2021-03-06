%% Photon count rate calculation from pre-calculated structural array containing position dependent

function [PCR] = PhotonCountRate( decayrates )
    
    %% Photon counts calculation
    PCR = struct();
    
    d = 1 ; % choose a fixed location for plotting of power dependent results. The closest points in d_BEM will be chosen according to decay rates results.
    
    %PCR.AbsorptionCrossSection = 10e-20 ; % m^2
    PCR.AbsorptionCrossSection = 5.7e-20 ; % m^2, ATTO647N 
    %PCR.AbsorptionCrossSection = 9.6e-20 ; % m^2, Cy5, CV
    PCR.collectionEfficiency = 0.02 ; % This value has been tune to match the experimental results
    PCR.h = 6.62e-34; % planck constant J * s
    PCR.lambda = decayrates.enei_field*1e-9 ; % m2
    PCR.c = 2.99792458e8 ; % m/s
    PCR.nu = PCR.c/PCR.lambda;
    %default lifetime of the dipole
    %lifetime = 1e-8; % seconds
    PCR.lifetime = 3.5e-9 ;% ATTO 647N
    %PCR.lifetime = 2.8e-10 ;% ATTO 647N, 13 times quenching, QY = 0.05
    %PCR.lifetime = 1e-9 ;% Cy5
    %PCR.lifetime = 0.1e-9; % Crystal violet, 100 ps ;
    %PCR.k_tot = 2.8e8; % no quenching
    %PCR.k_tot = 3.6e9; % 13 times quenching
    PCR.k_tot = 1/PCR.lifetime; % note! The non-enhanced total decay rate changes for intrinsic different quantum yields!
    
    PCR.k_r = 1.8e8; % decayrates.QY*PCR.k_tot; % no quenching
    %PCR.k_nr = 1.1e8; % PCR.k_tot - PCR.k_r; % no quenching
    %PCR.k_nr = 3.4e9; % 13 times quenching, QY = 0.05.
    PCR.k_nr = PCR.k_tot - PCR.k_r ;
    % set quenching by increasing the non-radiative rate. 
    %I_e = 1e-5; % Watt
    %I_inc = I_e ./( h*c/(enei_field*1e-9) * pi * ( 0.425 *(enei_field*1e-9)/NA )^2 ) *1e-18; % confocal photon density, photons nm-2 s-1
    
    PCR.I_inc = logspace(0,10,50) ; % incident power density, W m^-2, 1 - 1e10 W m^-2
    
    PCR.ISC_ratio = 300; % k_ISC/k_T. Wenger, Optics Express 2008
    triplet_enh = 1e4 ; % Use this term for hypothesized triplet decay enhancement. This term is effetively: (k_isc + k_T)/k_r^0.
    
    PCR.I_sat = PCR.h*PCR.nu*PCR.k_tot./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ;
    
    PCR.I_satenh = PCR.h*PCR.nu*(PCR.k_r.*decayrates.tot_average + PCR.k_r.*triplet_enh )./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ; % distance-dependent term from tot_average, this is a local field term. Far-field values needs to be devided by near field enhancement factor.
    % Note!
    % Modifed total decay rate is k_tot* = k_r * tot_average. Since tot_average
    % = k_tot*/(k_r)
    PCR.k_r_enh = PCR.k_r.*decayrates.rad_average; 
    PCR.k_nr_enh = PCR.k_r.*(decayrates.tot_average - decayrates.rad_average);
    PCR.k_tot_enh = PCR.k_r.*decayrates.tot_average; 
    
    %% Plot quantum yield and distance dependent saturation curves, at HIGH excitation power
    figure
    
    PCR.I_exc = PCR.I_inc(end) ; % w m^-2
    
    for i = 1 : length(decayrates.QY)
    %for i = 2   
        PCR.EnhCountRateDistQYdepend_HI(i,:)  = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*decayrates.QY(i).* decayrates.Q_avg(i,:)'./decayrates.QY(i).*PCR.I_exc.*PCR.I_satenh.*decayrates.ee./(PCR.I_exc.*decayrates.ee + PCR.I_satenh);% enhanced PCR count rate, distance dependent
        PCR.NonenhCountRateDistQYdepend_HI(i)  = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency....
            /(PCR.h*PCR.nu).*decayrates.QY(i).*PCR.I_exc.*PCR.I_sat./(PCR.I_exc + PCR.I_sat);% non-enhanced PCR count rate, non-distance dependent, fixed power
        plot( decayrates.d_BEM, PCR.EnhCountRateDistQYdepend_HI(i,:))
        hline(PCR.NonenhCountRateDistQYdepend_HI(i),':')
        
        hold on
    end
    xlabel('d (nm)')
    ylabel('Photon count rate (s^{-1})')
    set(gca,'YScale','log')
    title(['Laser power: ', num2str(PCR.I_exc,' %1.0e'),' W m^{-2}'])
    MinWhitSpace
    hold off
    
    saveas( gcf,'Distance QY dependent Photon_rates_HIGH_exc.fig' )
    saveas( gcf,'Distance QY dependent Photon_rates_HIGH_exc.png' )
    
    
    %% Plot quantum yield and distance dependent saturation curves, at LOW excitation power
    figure
    
    PCR.I_exc = PCR.I_sat ; % w m^-2
    
    for i = 1 : length(decayrates.QY)
     %for i = 2   
        PCR.EnhCountRateDistQYdepend_LO(i,:)  = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*decayrates.QY(i).* decayrates.Q_avg(i,:)'./decayrates.QY(i).*PCR.I_exc.*PCR.I_satenh.*decayrates.ee./(PCR.I_exc.*decayrates.ee + PCR.I_satenh);% enhanced PCR count rate, distance dependent
        PCR.NonenhCountRateDistQYdepend_LO(i)  = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency....
            /(PCR.h*PCR.nu).*decayrates.QY(i).*PCR.I_exc.*PCR.I_sat./(PCR.I_exc + PCR.I_sat);% non-enhanced PCR count rate, non-distance dependent, fixed power
        plot( decayrates.d_BEM, PCR.EnhCountRateDistQYdepend_LO(i,:))
        hline(PCR.NonenhCountRateDistQYdepend_LO(i),':')
        
        hold on
    end
    xlabel('d (nm)')
    ylabel('Photon count rate (s^{-1})')
    set(gca,'YScale','log')
    title(['Laser power: ', num2str(PCR.I_exc,' %1.0e'),' W m^{-2}'])
    MinWhitSpace
    hold off
    close
    saveas( gcf,'Distance QY dependent Photon_rates_LOW_exc.fig' )
    saveas( gcf,'Distance QY dependent Photon_rates_LOW_exc.png' )
    
    
    %% Plot power dependent saturation curves for all quantum yields at the chosen position (fixed at 3 nmor other)
    
    for i = 1 : length (decayrates.QY)
    %for i = 2
        %PCR.I_sat = PCR.h*PCR.nu*PCR.k_tot./PCR.AbsorptionCrossSection  ; % non-enhanced I_sat
        PCR.NonenhCountRatePowerQYdepend(i,:)  =  PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./(PCR.h*PCR.nu).*decayrates.QY(i).*PCR.I_inc.*PCR.I_sat./(PCR.I_inc + PCR.I_sat);% non-enhanced PCR count rate, non-distance dependent, power-dependent PCR
        
        %[~,idx] = max(FE(i,:));
        
        %idx = 23; % distance fixed at about 3 nm; idx = 23, distance at about 7 nm. 
        [~, dist_idx] = min(abs(decayrates.d_BEM - d)); % find the closest point to input distance and return index.
        
        PCR.I_incNF = PCR.I_inc.*decayrates.ee(dist_idx); % near field excitaiton intensity is calculated here for the enhanced photon rate calcualtion. This gives the far-field saturation intensity, which is reduced than non-enhanced fluorophores.
        PCR.EnhCountRatePowerQYdepend(i,:) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*decayrates.QY(i).*decayrates.Q_avg(i,dist_idx)./decayrates.QY(i).*PCR.I_incNF.*PCR.I_satenh(dist_idx)./(PCR.I_incNF + PCR.I_satenh(dist_idx)); % enhanced power dependent-PCR, note that local enhanced field is used here.
        
        %PCR.maxPCR(i) = interp1(PCR.I_inc,PCR.CountRatePD,2.*PCR.I_sat,'spline');
        PCR.maxPCR(i) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*decayrates.QY(i).* PCR.I_sat;  % non-enhanced maximal PCR
        PCR.maxPCRenh(i) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*decayrates.QY(i).*decayrates.Q_avg(i,dist_idx)./decayrates.QY(i).*PCR.I_satenh(dist_idx); % enhanced maximal PCR
        %PCR.maxPCRenh(i) = interp1(PCR.I_incNF,PCR.PCR_enh,2.*PCR.I_satenh(idx),'spline');
        PCR.maximumPCRenhancement(i) = PCR.maxPCRenh(i)./PCR.maxPCR(i);
        
        figure
        plot(PCR.I_inc, PCR.NonenhCountRatePowerQYdepend(i,:), PCR.I_inc, PCR.EnhCountRatePowerQYdepend(i,:))
        %annotation('textbox',[.2 .5 .3 .3],'String',{'I_sat = ',num2str(I_sat)},'FitBoxToText','on');
        hline(PCR.maxPCR(i),'-','Maximal PCR')
        hline(PCR.maxPCRenh(i),'-','Maximal PCR enhanced')
        vline(PCR.I_sat,':','Non-enhanced I_sat') % non_enhanced I_sat
        vline(PCR.I_satenh(dist_idx)./decayrates.ee(dist_idx),':','Enhanced I_sat')  % enhanced I_sat, converted back to far-field excitation intensity
        xlabel('Incident power (W m^{-2})')
        ylabel('Photon count rate (s^{-1})')
        title(['\Phi',num2str(decayrates.QY(i),'%1.2f'),' Maximal PCR enhancement: ', num2str(PCR.maxPCRenh(i)/PCR.maxPCR(i)),...
            ' I_{sat} enhancement:',num2str(PCR.I_satenh(dist_idx)./decayrates.ee(dist_idx)./PCR.I_sat)])
        
        MinWhitSpace
        saveas( gcf,['Power QY dependent Photon_rates QY ',num2str(decayrates.QY(i)),'.fig'] )
        saveas( gcf,['Power QY dependent Photon_rates QY ',num2str(decayrates.QY(i)),'.png'] )
        
    end
    %%
%     figure
%     loglog(decayrates.QY,PCR.maximumPCRenhancement,'ro-',decayrates.phi, decayrates.QY_mod, 'b-', decayrates.phi, decayrates.tot_mod, 'g-')
%     xlabel('\Phi_0')
%     ylabel('PCR^*_{max}/PCR_{max}')
%     legend('Maximal PCR enhancement','QY modification','lifetime modification')
%     MinWhitSpace
%     saveas( gcf,'QY dependent PCR max enhancement.fig' )
%     saveas( gcf,'QY dependent PCR max enhancement.png' )
% %% Plot the quantum yield depedent PCR_max and I_sat modifications. 
% %  This calculation is for showing only, and does not go into the result
% %  array.
%         %phi_0 = logspace( -2, 0, 1000 );
%         phi_0 = linspace(0.001, 1, 1000);
%         k_r = 1.8e8; % intrinsic radiative decay rate. 
%         k_tot = k_r./phi_0; % intrinsic total decay rate.
%     
%     for i = 1 : length(phi_0)
%         
%         phi(i,:) =  decayrates.rad_average./( decayrates.tot_average + (1 - phi_0(i) )./ phi_0(i) ) ;
%               
%         I_sat(i) = PCR.h*PCR.nu*k_tot(i)./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ;
%         
%         I_satenh(i,:) = PCR.h*PCR.nu*k_r.*decayrates.tot_average./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ; % distance-dependent term from tot_average, this is a local field term. Far-field values needs to be devided by near field enhancement factor.
%         
%         PCR_max(i,:) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
%             (PCR.h*PCR.nu).*phi_0(i).*phi(i,:)./phi_0(i).*I_satenh(i,:); % enhanced maximal PCR
%         PCR_maxnonenh(i) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
%             (PCR.h*PCR.nu).*phi_0(i).* I_sat(i);  % non-enhanced maximal PCR
%     end
%     
%     figure
%     plot(phi_0,PCR_max(:,dist_idx),phi_0,PCR_maxnonenh)
%     set(gca,'FontSize',10)
%     pbaspect([1.2, 1 , 1])
%     %gca.PlotBoxAspectRatio = [1 0.75 0.75];
%     %set(gcf,'units','centimeters','position',[10,10,14,10])
%     set(gca,'units','centimeters','position',[1,1,7,5])
%     
%     figure
%     loglog(phi_0,I_satenh(:,dist_idx), phi_0, I_sat, phi_0, I_satenh(:,dist_idx)./decayrates.ee(dist_idx) )
%     
%     figure;contourf(decayrates.d_BEM, diameter, I_satenhFFsort,1000, 'edgecolor', 'none')
%     xlabel('molecule-tip separation d (nm)')
%     ylabel('\phi_0')
%     legend('I_{sat} (W m^{-2})')
%     colormap jet
%     colorbar
%     set(gca,'colorscale', 'log', 'xscale', 'log')
%     set(gca,'FontSize',10)
%     pbaspect([1.2, 1 , 1])
%     %gca.PlotBoxAspectRatio = [1 0.75 0.75];
%     %set(gcf,'units','centimeters','position',[10,10,14,10])
%     set(gca,'units','centimeters','position',[3,3,7,5])
%     
%     % plot distance and phi_0 dependent phi
%     figure
%     [C,h] = contourf(decayrates.d_BEM,phi_0,phi);
%     xlabel('molecule-tip separation d (nm)')
%     ylabel('\phi_0')
%     %clabel(C,h,'FontSize',12)
%     %h.ShowText = 'on';
%     set(gca,'XScale','log','yscale','log')
%     colorbar
%     caxis([0.001 1])
%     colormap jet
%     
%     figure
%     imagesc(decayrates.d_BEM,phi_0,phi);
%     xlabel('molecule-tip separation d (nm)')
%     ylabel('\phi_0')
%     %clabel(C,h,'FontSize',12)
%     %h.ShowText = 'on';
%     %set(gca,'XScale','log','yscale','log')
%     colorbar
%     caxis([0.001 1])
%     colormap jet
%     
%     % plot distance and phi_0 dpendent PCR_max and I_sat
%     figure
%     [C,h] = contourf(decayrates.d_BEM, phi_0, PCR_max );
%     xlabel('molecule-tip separation d (nm)')
%     ylabel('\phi_0')
%     %clabel(C,h,'FontSize',10)
%     set(gca,'xscale','log','yscale','log','ColorScale','log')
%     %addlistener(h,'MarkedClean',@(a,b)ReFormatText(a))
%     colormap jet
%     colorbar
%     
%     figure
%     [C,h] = contourf(decayrates.d_BEM, phi_0, I_satenh);
%     xlabel('molecule-tip separation d (nm)')
%     ylabel('\phi_0')
%     %clabel(C,h,'FontSize',10)
%     set(gca,'xscale','log','yscale','log','ColorScale','log')
%     %addlistener(h,'MarkedClean',@(a,b)ReFormatText(a))
% 
%     colormap jet
%     colorbar
%% Plot radiative rate dependent PCR_max and I_sat. PCR_max and I_sat strongly depends on the intrinsic radiative rate. 
    
    phi_0 = 0.65; 
    
    k_r = logspace(4,9,100); % intrinsic radiative decay rate. 
    k_tot = k_r./phi_0; % intrinsic total decay rate.
    
    
    for i = 1 : length(k_r)
        
        phi(i,:) =  decayrates.rad_average./( decayrates.tot_average + (1 - phi_0)./ phi_0 ) ;
              
        I_sat(i) = PCR.h*PCR.nu*k_tot(i)./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ;
        
        I_satenh(i,:) = PCR.h*PCR.nu*k_r(i).*decayrates.tot_average./(PCR.AbsorptionCrossSection.*(1 + PCR.ISC_ratio))  ; % distance-dependent term from tot_average, this is a local field term. Far-field values needs to be devided by near field enhancement factor.
        
        PCR_max(i) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*phi_0.*phi(i,dist_idx)./phi_0.*I_satenh(i,dist_idx); % enhanced maximal PCR
        PCR_maxnonenh(i) = PCR.AbsorptionCrossSection.* PCR.collectionEfficiency./...
            (PCR.h*PCR.nu).*phi_0.* I_sat(i);  % non-enhanced maximal PCR
    end
    
    
    figure
    loglog(k_r,PCR_max,k_r,PCR_maxnonenh)
    set(gca,'FontSize',10)
    pbaspect([1.2, 1 , 1])
    %gca.PlotBoxAspectRatio = [1 0.75 0.75];
    %set(gcf,'units','centimeters','position',[10,10,14,10])
    set(gca,'units','centimeters','position',[1,1,7,5])
    
    figure
    loglog(k_r,I_satenh(:,dist_idx), k_r, I_sat, k_r, I_satenh(:,dist_idx)./decayrates.ee(dist_idx) )
    
%% Plot the ratio of lifetime reduction and near field enhancement as a distance dependent curve
    
    figure
    loglog(decayrates.d_BEM, decayrates.tot_average./decayrates.ee )
    hline(1,'-','No saturation intensity modifications')
    xlabel('d (nm)')
    ylabel('I^*_{sat}/{I_{sat}}')
    
    MinWhitSpace
    saveas( gcf,'DisDepend I_sat enhancemement.fig' )
    saveas( gcf,'DisDepend I_sat enhancemement.png' )
    
    figure
    loglog(decayrates.d_BEM, decayrates.tot_average(:),decayrates.d_BEM, decayrates.ee(:))
    xlabel('d (nm)')
    legend('k_{tot}^*/k_{tot}','|E|^2/|E_0|^2')
    MinWhitSpace
    saveas( gcf,'DisDepend k_tot enhancemement.fig' )
    saveas( gcf,'DisDepend k_tot enhancemement.png' )
    
    save PCR PCR
    
    %clearvars -except enei_field enei_dipole QY
    close all
    
end