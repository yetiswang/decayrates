%% Fluorescence enhancement of nanorods, results are calculated along the long axis
%  as a function of distance.
%  The script can calculate the overall fluorescence enhancement as a function of molecule
%  surface distance and dipole wavelength
%  Output: [FE] is the fluorescence enhancement factors as a function of distance and QY.
%                along row, distance-dependent total fluorescence enhancement
%                along column, QY in the input
function [PCR, FE, Lorentz, ee] = FE_GNR_1D( height, diameter, metal, enei_field, enei_dipole )

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

switch metal
    case 'Au'
        %  table of dielectric functions
        epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };
    case 'AgPalik'
        %  table of dielectric functions
        epstab = { epsconst( 1.33^2 ), epstable( 'silver_palik.dat' ) };
    case 'AgJC'
        %  table of dielectric functions
        epstab = { epsconst( 1.33^2 ), epstable( 'silver.dat' ) };
        
end


%  nanorod geometries
mesh = [ 41, 41, 41]; % n1 for the circumference of the rod, n2 for the polar angles of the rod caps, n3 for the cylinder-shaped middle part of the rod

QY = [ 0.01 0.02 0.05 0.1 0.2 0.5 1]; % series of QY for calculation


%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trirod ( diameter, height, mesh, 'triangles' );

%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);


%% make a new directory
directory = pwd;
switch metal
    case 'Au'
        ndir = ['H',num2str(height),'D',num2str(diameter),'_GNR_distance_dependent_FE','_Exc', num2str(enei_field),'_Dip_',num2str(enei_dipole)];
    case 'AgPalik'
        ndir = ['H',num2str(height),'D',num2str(diameter),'_AgNR_palik_distance_dependent_FE','_Exc', num2str(enei_field),'_Dip_',num2str(enei_dipole)];
    case 'AgJC'
        ndir = ['H',num2str(height),'D',num2str(diameter),'_AgNR_JC_distance_dependent_FE','_Exc', num2str(enei_field),'_Dip_',num2str(enei_dipole)];
end
mkdir(ndir)
cd(ndir)

%% Calculate and save the scattering spectrum for reference
[ ~, ~, Lorentz ] = spect_GNR_BEM( epstab, height, diameter, linspace(500,1000,50 ));

%%  set wavelength of planewave and dipole oscillator
if nargin == 3
    enei_field = 1248./Lorentz(3) ;  % 637
    enei_dipole = enei_field + 25 ; % 670
else
end
vline(enei_dipole, 'r-','\lambda_{dip}');
vline(enei_field,  'b-','\lambda_{exc}');
saveas (gcf, [ndir,'.fig'], 'fig')
saveas (gcf, [ndir,'.png'], 'png')
%% Dipole positions

%1D positions of dipole
%choose the minimum distance by setting the values of x vector
x = reshape( logspace( log10(0.51 * height), log10( 100 + 0.5 * height ), 200 ), [], 1 );
%x = reshape( linspace( 0.51, (50 + 0.5 * height )/height, 100 ) * height, [], 1 );

% compoint
pt = compoint( p, [ x, 0 .* x, 0 .* x ], 'mindist' , 1e-10  );

dir_dip = [ 1, 0, 0 ; 0 , 1 , 0 ; 0, 0, 1];

% dipole excitation , x and z direction
dip = dipole( pt, dir_dip, op );

% plot dipole orientation and nanorod
figure
plot(p)
axis on
hold on
plot3(pt.pos(:,1), pt.pos(:,2),pt.pos(:,3) ,'r.')
xlabel('x(nm)')
ylabel('z(nm)')
hold off
view([0 0])
saveas(gcf,'dipole-particle.fig')
saveas(gcf,'dipole-particle.png')
%%  BEM simulation for decay rates
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei_dipole );
%  total and radiative decay rate
[ tot, rad, ~] = dip.decayrate( sig );
%% decay rate plot for intrinsic quantum yield of 1%. These plots are generated as examples and previews.
d_BEM = x - height/2;

figure
semilogy( x , tot, '-'  );  hold on;
semilogy( x, rad, 'o-' );
xlim( [ min( x ), max( x ) ] );
title( 'Total and radiaitve decay rate for dipole oriented along x and z' );
legend( 'tot_x BEM','tot_y BEM','tot_z BEM','rad_x BEM','rad_y BEM', 'rad_z BEM' )

xlabel( 'Position (nm)' );
ylabel( 'Decay rate' );
MinWhitSpace
saveas(gcf,'Decay rates.fig')
saveas(gcf,'Decay rates.png')

Qx = rad(:, 1)./( tot(:, 1) + (1 - 0.01 ) / 0.01);
Qy = rad(:, 2)./( tot(:, 2) + (1 - 0.01 ) / 0.01);
Qz = rad(:, 3)./( tot(:, 3) + (1 - 0.01 ) / 0.01);

rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;

Q_exam =  rad_average./( tot_average + (1 - 0.01 ) / 0.01 ) ;

figure
plot(d_BEM, Qx , d_BEM, Qy, d_BEM, Qz, d_BEM, Q_exam  )
legend('QY_x BEM', 'QY_y BEM', 'QY_z BEM','QY_average BEM')
xlabel('Distance to surface (nm)')
ylabel(['Quantum yield (intrinsic QY = ', num2str(0.01), ' )'])
title(['Wavelength of dipole ', num2str(enei_dipole),' nm'])
xlim( [ 0 , max(d_BEM) ]  )
MinWhitSpace
saveas(gcf,'QY_axis.fig')
saveas(gcf,'QY_axis.png')

%% BEM solver for near field intensity

%  TM mode, excitation from above
dir = [ 0, 0, -1];
pol = [ 1, 0, 0 ];

%  initialize BEM solver
bem = bemsolver( p, op );
%  initialize plane wave excitation
exc = planewave( pol, dir, op );
%  solve BEM equation
sig = bem \ exc( p, enei_field );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  close waitbar
multiWaitbar( 'CloseAll' );

% %  set up Green function object between PT and P
% %  use the pt object calculated from above dipole positions
% g = greenfunction( pt, p, op );
% %  compute electric field
% f = field( g, sig );



% ee = sqrt( dot (f.e , f.e, 2 ) ) ;
% ee = ee.^2 ;
% ee_norm = vecnorm( f.e ) ;
% ee_normsquared = ee_norm.^2;

%%  computation of electric field
%  object for electric field
%    MINDIST controls the minimal distance of the field points to the
%    particle boundary, MESHFIELD must receive the OP structure which also
%    stores the table of precomputed reflected Green functions
emesh = meshfield( p, x, 0 .* x, 0 .* x, op, 'mindist', 0.2, 'nmax', 2000 , 'waibar', 1 );
%  induced and incoming electric field
e = emesh( sig ) + emesh( exc.field( emesh.pt, enei_field ) );
%  norm of electric field
enorm = vecnorm(e);

ee = enorm.^2;

% plot near field
% get electric field
ex = sqrt( dot (e(: ,1) , e(: ,1), 	2 ) );
ey = sqrt( dot (e(: ,2) , e(: ,2), 	2 ) );
ez = sqrt( dot (e(: ,3) , e(: ,3), 	2 ) );


%e_ind = emesh(sig);
%e_ext = emesh(exc.field( emesh.pt, enei_field ));
%%
% plot enhanced field in averaged orietation
figure
loglog(d_BEM ,ee(:))
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity')
MinWhitSpace
saveas(gcf, 'nearfield_axis.fig')
saveas(gcf, 'nearfield_axis.png')

% % plot electric field vector
% figure
% coneplot( pt.pos, f.e )
% axis on
% grid on
% hold on;
% plot(p)
% saveas( gcf, 'Electric field vector.fig' )
% saveas( gcf, 'Electric field vector.png' )

% plot enhanced field in every component
figure
loglog(d_BEM, ex.^2,d_BEM, ey.^2,d_BEM, ez.^2);
legend('E_x^2/E_0^2','E_y^2/E_0^2','E_z^2/E_0^2')
xlabel('d (nm)')
ylabel('Near field enhancement')
MinWhitSpace
saveas( gcf, 'nearfield_xyz.fig' )
saveas( gcf, 'nearfield_xyz.png' )

%% Combine QY and NF, calculate curve for different positions
rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;

for i = 1 : length(QY)
    
    Q_avg(i,:) =  rad_average./( tot_average + (1 - QY(i) ) / QY(i) ) ;
    FE(i,:) = ee(:)'.*Q_avg(i,:)./QY(i) ;
    
    
    figure
    semilogy(d_BEM, FE(i,:), 'r-o','LineWidth',1)
    legend('\xi')
    xlabel('d (nm)')
    ylabel('Fluorescence enhancement')
    str = ['\lambda_{dipole} = ', num2str(enei_dipole),' nm \newline \lambda_{exc} = ', ...
        num2str(enei_field), ' nm'];
    annotation('TextBox',[0.5 0.3 0.5 0.6],'String',str,'FitBoxToText','on','LineStyle','none','FontSize',14);
    MinWhitSpace
    saveas(gcf,[num2str(QY(i)),'_fluorescence_enhancement_on_axis.fig'])
    saveas(gcf,[num2str(QY(i)),'_fluorescence_enhancement_on_axis.png'])
    
end
%% save data
save d_BEM d_BEM
save gamma_rad_BEM rad
save gamma_tot_BEM tot
save nearfield ee
save nearfieldx ex
save nearfieldy ey
save nearfieldz ez
save FE_QY FE
save field e
% save FEx FEx
% save FEy FEy
% save FEz FEz
save FE FE
save Qavg Q_avg

%% Plot QY vs Max. FE

phi = linspace( 0, 1, 1000 );

for i = 1 : length(phi)
    
    Q_line =  rad_average./( tot_average + (1 - phi(i) )./ phi(i) ) ;
    
    factor = Q_line./phi(i).*ee;
    
    %[FE_allphi(i),idx] = max(factor);
    idx = 11; % fixed distance 3 nm. 
    tot_mod(i) = tot_average(idx);
    
    QY_mod(i) = Q_line(idx)./phi(i);
    
end
%% 

figure
semilogx( phi, FE_allphi )
xlabel('QY')
ylabel('\xi_{EF}')
hline(max(ee),'r-','max. exc. enhancement')
MinWhitSpace
saveas(gcf,'MaxFEvsphi.fig')
saveas(gcf,'MaxFEvsphi.png')

save FE_allphi FE_allphi
save phi phi


%% Photon counts calculation
PCR = struct();

AbsorptionCrossSection = 10e-20 ; % m^2
AbsorptionCrossSection = 1e-20 ; % m^2
collectionEfficiency = 0.08 ;
h = 6.62e-34; % planck constant J * s
lambda = enei_field*1e-9 ; % m2
c = 2.99792458e8 ; % m/s
nu = c/lambda;
% default lifetime of the dipole
lifetime = 1e-8; % seconds
lifetime = 3.5e-9 ;% 
k_tot = 1/lifetime;
%I_e = 1e-5; % Watt
%I_inc = I_e ./( h*c/(enei_field*1e-9) * pi * ( 0.425 *(enei_field*1e-9)/NA )^2 ) *1e-18; % confocal photon density, photons nm-2 s-1

PCR.I_inc = logspace(0,7.5,50) ; % incident power density, W m^-2, 1 - 1e10 W m^-2
PCR.I_satenh = h*nu*k_tot.*tot_average./AbsorptionCrossSection  ; % distance-dependent term from tot_average, this is a local field term. Far-field values needs to be devided by near field enhancement factor.

%% Plot quantum yield and distance dependent saturation curves
figure

PCR.I_exc = 1e8 ; % w m^-2
PCR.I_sat = h*nu*k_tot./AbsorptionCrossSection  ;
for i = 1 : length(QY)
    %Q_avg =  rad_average./( tot_average + (1 - QY(i) )./ QY(i) ) ;
    PCR.EnhCountRateDistQYdepend(i,:)  = AbsorptionCrossSection.* collectionEfficiency./...
        (h*nu).*QY(i).* Q_avg(i,:)'./QY(i).*PCR.I_exc.*PCR.I_satenh.*ee./(PCR.I_exc.*ee + PCR.I_satenh);% enhanced PCR count rate, distance dependent
    PCR.NonenhCountRateDistQYdepend(i)  = AbsorptionCrossSection.* collectionEfficiency....
        /(h*nu).*QY(i).*PCR.I_exc.*PCR.I_sat./(PCR.I_exc + PCR.I_sat);% non-enhanced PCR count rate, non-distance dependent, fixed power
    plot( d_BEM, PCR.EnhCountRateDistQYdepend(i,:))
    hline(PCR.NonenhCountRateDistQYdepend(i),':')
    
    hold on
end
xlabel('d (nm)')
ylabel('Photon count rate (s^{-1})')
set(gca,'YScale','log')
title(['Laser power: ', num2str(PCR.I_exc,' %1.0e'),' W m^{-2}'])
MinWhitSpace
hold off
saveas( gcf,'Distance QY dependent Photon_rates.fig' )
saveas( gcf,'Distance QY dependent Photon_rates.png' )

%% Plot quantum yield and power dependent saturation curves at the optimal position (fixed at 3 nm)

for i = 1 : length (QY)
    PCR.I_sat = h*nu*k_tot./AbsorptionCrossSection  ; % non-enhanced I_sat
    PCR.NonenhCountRatePowerQYdepend(i,:)  =  AbsorptionCrossSection.* collectionEfficiency./(h*nu).*QY(i).*PCR.I_inc.*PCR.I_sat./(PCR.I_inc + PCR.I_sat);% non-enhanced PCR count rate, non-distance dependent, power-dependent PCR
    
    %[~,idx] = max(FE(i,:));
    
    idx = 11; % distance fixed at about 3 nm;
    
    PCR.I_incNF = PCR.I_inc.*ee(idx); % near field excitaiton intensity is calculated here for the enhanced photon rate calcualtion. This gives the far-field saturation intensity, which is reduced than non-enhanced fluorophores.
    PCR.EnhCountRatePowerQYdepend(i,:) = AbsorptionCrossSection.* collectionEfficiency./...
        (h*nu).*QY(i).*Q_avg(i,idx)./QY(i).*PCR.I_incNF.*PCR.I_satenh(idx)./(PCR.I_incNF + PCR.I_satenh(idx)); % enhanced power dependent-PCR, note that local enhanced field is used here.
    
    %PCR.maxPCR(i) = interp1(PCR.I_inc,PCR.CountRatePD,2.*PCR.I_sat,'spline');
    PCR.maxPCR(i) = AbsorptionCrossSection.* collectionEfficiency./...
        (h*nu).*QY(i).* PCR.I_sat;  % non-enhanced maximal PCR
    PCR.maxPCRenh(i) = AbsorptionCrossSection.* collectionEfficiency./...
        (h*nu).*QY(i).*Q_avg(i,idx)./QY(i).*PCR.I_satenh(idx); % enhanced maximal PCR
    %PCR.maxPCRenh(i) = interp1(PCR.I_incNF,PCR.PCR_enh,2.*PCR.I_satenh(idx),'spline'); 
    PCR.maximumPCRenhancement(i) = PCR.maxPCRenh(i)./PCR.maxPCR(i);
    
    figure
    plot(PCR.I_inc, PCR.NonenhCountRatePowerQYdepend(i,:), PCR.I_inc, PCR.EnhCountRatePowerQYdepend(i,:))
    %annotation('textbox',[.2 .5 .3 .3],'String',{'I_sat = ',num2str(I_sat)},'FitBoxToText','on');
    hline(PCR.maxPCR(i),'-','Maximal PCR')
    hline(PCR.maxPCRenh(i),'-','Maximal PCR enhanced')
    vline(PCR.I_sat,':','Non-enhanced I_sat') % non_enhanced I_sat
    vline(PCR.I_satenh(idx)./ee(idx),':','Enhanced I_sat')  % enhanced I_sat, converted back to far-field excitation intensity
    xlabel('Incident power (W m^{-2})')
    ylabel('Photon count rate (s^{-1})')
    title(['\Phi',num2str(QY(i),'%1.2f'),' Maximal PCR enhancement: ', num2str(PCR.maxPCRenh(i)/PCR.maxPCR(i)),...
        ' I_{sat} enhancement:',num2str(PCR.I_satenh(idx)./ee(idx)./PCR.I_sat)])
    
    MinWhitSpace
    saveas( gcf,['Power QY dependent Photon_rates QY ',num2str(QY(i)),'.fig'] )
    saveas( gcf,['Power QY dependent Photon_rates QY ',num2str(QY(i)),'.png'] )
    
end
figure
loglog(QY,PCR.maximumPCRenhancement,'ro-',phi, QY_mod, 'b-', phi, tot_mod, 'g-')
xlabel('\Phi_0')
ylabel('PCR^*_{max}/PCR_{max}')
legend('Maximal PCR enhancement','QY modification','lifetime modification')
MinWhitSpace
saveas( gcf,'QY dependent PCR max enhancement.fig' )
saveas( gcf,'QY dependent PCR max enhancement.png' )

%% Plot the ratio of lifetime reduction and near field enhancement as a distance dependent curve

figure
loglog(d_BEM, tot_average./ee )
hline(1,'-','No saturation intensity modifications')
xlabel('d (nm)')
ylabel('I^*_{sat}/{I_{sat}}')

MinWhitSpace
saveas( gcf,'DisDepend I_sat enhancemement.fig' )
saveas( gcf,'DisDepend I_sat enhancemement.png' )

figure
loglog(d_BEM, tot_average(:),d_BEM, ee(:))
xlabel('d (nm)')
legend('k_{tot}^*/k_{tot}','|E|^2/|E_0|^2')
MinWhitSpace
saveas( gcf,'DisDepend k_tot enhancemement.fig' )
saveas( gcf,'DisDepend k_tot enhancemement.png' )

save PCR PCR

clearvars -except enei_field enei_dipole QY 
close all
end
