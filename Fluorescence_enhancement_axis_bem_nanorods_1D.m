%% Fluorescence enhancement of nanorods, results are calculated along the long axis
%  as a function of distance.
%  The script can calculate the overall fluorescence enhancement as a function of molecule
%  surface distance and dipole wavelength

clear all
close all
clc

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  table of dielectric functions
epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };

%  nanorod geometries
mesh = [ 41, 41, 41]; % n1 for the circumference of the rod, n2 for the polar angles of the rod caps, n3 for the cylinder-shaped middle part of the rod
height =  75  ;
diameter =  30  ;

Q_int = 0.01 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trirod ( diameter, height, mesh, 'triangles' );

%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);

%% Calculate and save the scattering spectrum for reference
[ sca, fit, Lorentz ] = spect_GNR_BEM( epstab, height, diameter, linspace(400,950,50 ));

%%  set wavelength of planewave and dipole oscillator

%enei_field = 1248./Lorentz(3) ;  % 637
enei_field = 610 ;
%enei_field = 532 ;
enei_dipole = enei_field + 25 ; % 670
enei_dipole = 650 ;
%% make a new directory
directory = pwd;
ndir = ['H',num2str(height),'D',num2str(diameter),'_GNR_distance_dependent_FE','_Exc',num2str(enei_field),'_Dip',num2str(enei_dipole)];
mkdir(ndir)
cd(ndir)

vline(enei_dipole, 'r-','\lambda_{dip}');
vline(enei_field,  'b-','\lambda_{exc}');
saveas (gcf, [ndir,'.fig'], 'fig')
saveas (gcf, [ndir,'.png'], 'png')
%% Dipole positions

%1D positions of dipole
%choose the minimum distance by setting the values of x vector
x = reshape( logspace( log10(0.51 * height), log10( 50 + 0.5 * height ), 100 ), [], 1 );
%z = reshape( logspace( log10(0.51 * diameter), log10( 50 + 0.5 * diameter ), 100 ), [], 1 );
%x = reshape( linspace( 0.51, (50 + 0.5 * height )/height, 100 ) * height, [], 1 );

% compoint
pt = compoint( p, [ x, 0 .* x, 0 .* x ], 'mindist' , 1e-10  );

%pt = compoint( p, [ 0.* z, 0 .* z, z ], 'mindist' , 1e-10  );

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
[ tot, rad, rad0 ] = dip.decayrate( sig );
%%  final plot
d_BEM = x - height/2;
% d_BEM = z - diameter/2;

figure
semilogy( x , tot, '-'  );  hold on;
semilogy( x, rad, 'o-' );
xlim( [ min( x ), max( x ) ] );
% semilogy( z , tot, '-'  );  hold on;
% semilogy( z, rad, 'o-' );
% xlim( [ min( z ), max( z ) ] );
title( 'Total and radiaitve decay rate for dipole oriented along x and z' );
legend( 'tot_x BEM','tot_y BEM','tot_z BEM','rad_x BEM','rad_y BEM', 'rad_z BEM' )
xlabel( 'Position (nm)' );
ylabel( 'Decay rate' );

saveas(gcf,'Decay rates.fig')
saveas(gcf,'Decay rates.png')

Qx = rad(:, 1)./( tot(:, 1) + (1 - Q_int ) / Q_int);
Qy = rad(:, 2)./( tot(:, 2) + (1 - Q_int ) / Q_int);
Qz = rad(:, 3)./( tot(:, 3) + (1 - Q_int ) / Q_int);

rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;

Q_avg =  rad_average./( tot_average + (1 - Q_int ) / Q_int ) ;

figure
plot(d_BEM, Qx , d_BEM, Qy, d_BEM, Qz, d_BEM, Q_avg  )
legend('QY_x BEM', 'QY_y BEM', 'QY_z BEM','QY_average BEM')
xlabel('Distance to surface (nm)')
ylabel(['Quantum yield (intrinsic QY = ', num2str(Q_int), ' )'])
title(['Wavelength of dipole ', num2str(enei_dipole),' nm'])
xlim( [ 0 , max(d_BEM) ]  )
saveas(gcf,'QY_axis.fig')
saveas(gcf,'QY_axis.png')
%% Calculate near field intensity on the dipole positions


%%   BEM solver for near field intensity

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

%  set up Green function object between PT and P
%  use the pt object calculated from above dipole positions
g = greenfunction( pt, p, op );
%  compute electric field
f = field( g, sig );

%% plot near field
%  get electric field
ex = sqrt( dot (f.e(: ,1) , f.e(: ,1), 	2 ) );
ey = sqrt( dot (f.e(: ,2) , f.e(: ,2), 	2 ) );
ez = sqrt( dot (f.e(: ,3) , f.e(: ,3), 	2 ) );

ee = sqrt( dot (f.e , f.e, 2 ) ) ;
ee = ee.^2 ;

% plot enhanced field in averaged orietation
figure
plot(d_BEM ,ee(:))
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity')
saveas(gcf, 'nearfield_axis.fig')
saveas(gcf, 'nearfield_axis.png')

% plot electric field vector
figure
coneplot( pt.pos, f.e )
axis on
grid on
hold on;
plot(p)
saveas( gcf, 'Electric field vector.fig' )
saveas( gcf, 'Electric field vector.png' )

% plot enhanced field in every component
figure
semilogy(d_BEM, ex.^2,d_BEM, ey.^2,d_BEM, ez.^2);
legend('E_x^2/E_0^2','E_y^2/E_0^2','E_z^2/E_0^2')
xlabel('d (nm)')
ylabel('Near field enhancement')
saveas( gcf, 'nearfield_xyz.fig' )
saveas( gcf, 'nearfield_xyz.png' )

%% Combine QY and NF, calculate curve for different positions

QY = [ 0.01, 0.05, 0.25, 0.4, 0.75 ];

for i = 1 : length(QY)
    
    rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
    tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;
    Q_avg =  rad_average./( tot_average + (1 - QY(i) ) / QY(i) ) ;
    
    FE(i,:) = ee.*Q_avg./QY(i) ;
    
    figure
    plot(d_BEM, FE(i,:), 'r-o','LineWidth',1)
    legend('\xi')
    xlabel('d (nm)')
    ylabel('Fluorescence enhancement')
    str = ['\lambda_{dipole} = ', num2str(enei_dipole),' nm \newline \lambda_{exc} = ', ...
        num2str(enei_field), ' nm'];
    h= annotation('TextBox',[0.5 0.3 0.5 0.6],'String',str,'FitBoxToText','on','LineStyle','none','FontSize',14);
    saveas(gcf,[num2str(QY(i)),'_fluorescence_enhancement_on_axis.fig'])
    saveas(gcf,[num2str(QY(i)),'_fluorescence_enhancement_on_axis.png'])
    
    mkdir(['QY',num2str(QY(i))])
    cd(['QY',num2str(QY(i))])
    save Qavg Q_avg
    
    cd(directory)
    cd(ndir)
    
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
% save FEx FEx
% save FEy FEy
% save FEz FEz
save FE FE


%% Plot QY vs Max. FE

QY = linspace( 0, 1, 1000 );

for i = 1 : length(QY)
    
    Q_avg =  rad_average./( tot_average + (1 - QY(i) )./ QY(i) ) ;
    
    factor = Q_avg./QY(i).*ee;
    
    FE_allQY(i) = max(factor);
    
end
figure
semilogx(QY,FE_allQY)
xlabel('QY')
ylabel('\xi_{EF}')
hline(max(ee),'r-','max. exc. enhancement')
saveas(gcf,'FEvsQY.fig')
saveas(gcf,'FEvsQY.png')

save FE_allQY FE_allQY
save QY QY


plot( d_BEM, ee, d_BEM, Q_avg./Q_int, d_BEM,FE(1,:)  )
xlabel('d (nm)')
legend( '\xi_{exc}(|E|^2)','\xi_{em}','\xi_{FE}' )
set(gca,'YScale','log')
pbaspect([1,1,1])
saveas(gcf,'ALL_PLOT_QY0.01.fig')
saveas(gcf,'ALL_PLOT_QY0.01.png')


%% Photon counts calculation
QY = [ 0.01 0.05 0.1 0.25 0.75 0.9 ];
rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;
tau_bm = 1e-8; % s
n_w = 1.33;
NA = 1.2;
AbsorptionCrossSection = 0.01 ; % nm^-2
collectionEfficiency = 0.06 ;
h = 6.62e-34; % planck constant J * s
c = 3e8 ; % m/s
% default lifetime of the dipole
%lifetime = 1e-8; % seconds
lifetime = 6e-10; % Rachel's Atto 740
k_tot = 1/lifetime;
I_e = 1e-5; % Watt
I_inc = I_e ./( h*c/(enei_field*1e-9) * pi * ( 0.425 *(enei_field*1e-9)/NA )^2 ) *1e-18; % confocal photon density, photons nm-2 s-1
kappa = 2.33 * n_w / NA ;
radius = (enei_field/( 2 * NA)) /2.35; % lateral resolution. FWHM = 2.35 * sigma, if the diffraction limited spot is the FWHM, then the sigma is FWHM/2.35
volume = pi^(3/2) * kappa * ( radius*2.35*0.85)^3 ; % diffraction limited spot
centerGs = [ 0 0 0 ];
params = [0 1 centerGs radius*1e-9 radius*1e-9 radius*kappa*1e-9];
I_saturation = k_tot.*tot_average./AbsorptionCrossSection  ;

I_saturation_bg = k_tot./AbsorptionCrossSection  ;

figure
for i = 1 : length(QY)
    Q_avg =  rad_average./( tot_average + (1 - QY(i) )./ QY(i) ) ;
    I_ph(i,:)  = QY(i).* Q_avg./QY(i).*I_inc.*I_saturation.*ee./(I_inc.*ee + I_saturation).* AbsorptionCrossSection.* collectionEfficiency ;
    I_ph_bg(i,:) = QY(i).* I_inc.*I_saturation_bg./(I_inc + I_saturation_bg).* AbsorptionCrossSection * collectionEfficiency ;

    plot( d_BEM, I_ph(i,:) )
    hold on
    hline(I_ph_bg(i),':')
    
end
legend( 'QY 0.01','QY 0.05','QY 0.1','QY 0.25','QY 0.75' )
set(gca,'YScale','log')
title(['Laser power', num2str(I_e)])

hold off
saveas( gcf,'Photon_rates.fig' )
saveas( gcf,'Photon_rates.png' )
