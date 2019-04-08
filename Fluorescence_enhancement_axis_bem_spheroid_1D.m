%% Fluorescence enhancement of spheroids, results are calculated along the long axis
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
epstab = { epsconst( 1.5^2 ), epstable( 'silver_palik.dat' ) };

%  nanorod geometries
mesh = 2 ; 
height =  110   ;
diameter =  44  ;

Q_int = 0.01 ; % intrinsic quantum efficiency

p = SphrdBEM ( height, diameter, mesh  );

%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%% Calculate and save the scattering spectrum for reference
[ sca, fit, Lorentz ] = spect_Spheroid_BEM( epstab, height, diameter, linspace(600,1200,20 ));

%%  set wavelength of planewave and dipole oscillator

enei_field = 1248./Lorentz(3) ;  % 637
enei_dipole = enei_field + 25 ; % 670

%% Choose dependencies
wavelength_depedent = 0;
distance_dependent = 1;

%% make a new directory
directory = pwd;
ndir = ['H',num2str(height),'D',num2str(diameter),'_AgSpheroid_distance_dependent_decayrates','_Exc',num2str(enei_field,'%.1f'),'_Dip',num2str(enei_dipole,'%.1f')];
mkdir(ndir)
cd(ndir)


vline(enei_dipole, 'r-','\lambda_{dip}');
vline(enei_field,  'b-','\lambda_{exc}');
saveas (gcf, [ndir,'.fig'], 'fig')
saveas (gcf, [ndir,'.png'], 'png')
%% 
if distance_dependent == 1

%1D positions of dipole
%choose the minimum distance by setting the values of x vector
x = reshape( logspace( log10(0.51 * height), log10( 50 + 0.5 * height ), 100 ), [], 1 );
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
[ tot, rad, rad0 ] = dip.decayrate( sig );
%%  final plot
d_BEM = x - height/2;

figure
semilogy( x , tot, '-'  );  hold on;
semilogy( x, rad, 'o-' );
xlim( [ min( x ), max( x ) ] );
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

%% Combine QY and NF

% % QY enhancement calculated from QY of each component
% QYEx = Qx./Q_int ;
% QYEy = Qy./Q_int ;
% QYEz = Qz./Q_int ;
% % Fluorescence enhancement
% FEx = ex.^2.*QYEx ; 
% FEy = ey.^2.*QYEy ;
% FEz = ez.^2.*QYEz ;
FE = ee.*Q_avg./Q_int ;

% plot FE as a function of distance
% figure
% semilogy(d_BEM, FEx, d_BEM, FEy, d_BEM, FEz, d_BEM, FEavg, 'r-o','LineWidth',1)
% legend( 'F_x','F,y','F_z','F_{average}' )
% xlabel('d (nm)')
% ylabel('Fluorescence enhancement')
% str = ['\lambda_{dipole} = ', num2str(enei_dipole),' nm \newline \lambda_{exc} = ', ... 
%     num2str(enei_field), ' nm'];
% h= annotation('TextBox',[0.5 0.3 0.5 0.6],'String',str,'FitBoxToText','on','LineStyle','none','FontSize',14);
figure
plot(d_BEM, FE, 'r-o','LineWidth',1)
legend('\xi')
xlabel('d (nm)')
ylabel('Fluorescence enhancement')
str = ['\lambda_{dipole} = ', num2str(enei_dipole),' nm \newline \lambda_{exc} = ', ... 
    num2str(enei_field), ' nm'];
h= annotation('TextBox',[0.5 0.3 0.5 0.6],'String',str,'FitBoxToText','on','LineStyle','none','FontSize',14);
saveas(gcf,'fluorescence_enhancement_on_axis.fig')
saveas(gcf,'fluorescence_enhancement_on_axis.png')

%% save data
save d_BEM d_BEM
save gamma_rad_BEM rad
save gamma_tot_BEM tot
save nearfield ee
save nearfieldx ex
save nearfieldy ey
save nearfieldz ez
save Qx Qx
save Qy Qy
save Qz Qz
save Qavg Q_avg
% save FEx FEx
% save FEy FEy
% save FEz FEz
save FE FE

else
end 

if wavelength_depedent == 1
 
%%  dipole oscillator
    enei_dipole = linspace( 400, 900, 50 );
	enei_field = 637 ;
else
end