%% Fluorescence enhancement of nanospheres, results are calculated along the long axis
%  as a function of distance.
%  The result of this script is used for illustration

function FE_GNS_1D( diameter, metal )

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

Q_int = 0.01 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 51 ),  ...
    pi * linspace( 0, 1, 51 ) .^ 2, diameter, 'triangles' );

% initialize particle
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );



%% make a new directory
directory = pwd;
switch metal
    case 'Au'
        ndir = ['D',num2str(diameter),'_GNS_distance_dependent_FE'];
    case 'AgPalik'
        ndir = ['D',num2str(diameter),'_AgNS_palik_distance_dependent_FE'];
    case 'AgJC'
        ndir = ['D',num2str(diameter),'_AgNS_JC_distance_dependent_FE'];
end
mkdir(ndir)
cd(ndir)
%% Calculate and save the scattering spectrum for reference
[ ~, ~, Lorentz ] = spect_GNS_BEM( epstab, diameter, linspace(500,700,20 ));

%% Set wavelengths of dipole and excitation
enei_field = 1248./Lorentz(3) ;  % 637
enei_dipole = enei_field + 25 ; % 670
vline(enei_dipole, 'r-','\lambda_{dip}');
vline(enei_field,  'b-','\lambda_{exc}');
saveas (gcf, [ndir,'.fig'], 'fig')
saveas (gcf, [ndir,'.png'], 'png')
%%  dipole oscillator

%1D positions of dipole
%choose the minimum distance by setting the values of x vector
% x = reshape( linspace( 0.515, (50 + 0.5 * diameter )/diameter, 50 ) * diameter, [], 1 );
z = reshape( logspace( log10(0.51 * diameter), log10((50 + 0.5 * diameter )), 100 ) , [], 1 );

% compoint
pt = compoint( p, [ 0 * z, 0 * z, z ], 'mindist' , 1e-10  );

% dipole excitation , x y  and x direction
dip = dipole( pt, [ 1, 0, 0 ; 0 , 1 , 0 ; 0, 0 , 1], op );

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
d_BEM = pt.pos(:,3) - diameter/2;

figure
semilogy( z , tot, '-'  );  hold on;
semilogy( z, rad, 'o-' );
xlim( [ min( z ), max( z ) ] );
title( 'Total and radiaitve decay rate for dipole oriented along x, y and z' );
legend( 'tot_x BEM','tot_y BEM','tot_z BEM','rad_x BEM','rad_y BEM','rad_z BEM' )
xlabel( 'Position (nm)' );
ylabel( 'Decay rate' );

saveas(gcf,'decayrates.fig')
saveas(gcf,'decayrates.png')

Qx = rad(:, 1)./( tot(:, 1) + (1 - Q_int ) / Q_int);
Qy = rad(:, 2)./( tot(:, 2) + (1 - Q_int ) / Q_int);
Qz = rad(:, 3)./( tot(:, 3) + (1 - Q_int ) / Q_int);

rad_average = ( rad(:, 1) + rad(:, 2) + rad(:, 3) )./3 ;
tot_average = ( tot(:, 1) + tot(:, 2) + tot(:, 3) )./3 ;

Q_avg =  rad_average./( tot_average + (1 - Q_int ) / Q_int ) ;

figure
plot(d_BEM, Qx , d_BEM, Qy , d_BEM, Qz, d_BEM, Q_avg )
legend('QY_x BEM', 'QY_y BEM','QY_z BEM','QY_average BEM')
xlabel('Distance to surface (nm)')
ylabel(['Quantum yield (intrinsic QY = ', num2str(Q_int), ' )'])
title(['Wavelength of dipole ', num2str(enei_dipole),' nm'])
xlim( [ 0 , max(d_BEM) ]  )
hline(Q_int, 'r:','Intrisic QY')
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


pt_field = compoint( p, [ z, 0 * z, 0 * z ], 'mindist' , 1e-10  );

%  set up Green function object between PT and P
%  use the pt object calculated from above dipole positions
g = greenfunction( pt_field, p, op );
%  compute electric field
f = field( g, sig );

%% plot near field
%  plot electric field
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

%% plot enhanced field in every component
figure
semilogy(d_BEM, ex.^2,d_BEM, ey.^2,d_BEM, ez.^2);
legend('E_x^2/E_0^2','E_y^2/E_0^2','E_z^2/E_0^2')
xlabel('d (nm)')
ylabel('Near field enhancement')
saveas( gcf, 'nearfield_xyz.fig' )
saveas( gcf, 'nearfield_xyz.png' )


%% Combine QY and NF, calculate curve for different positions

QY = [ 0.01, 0.05, 0.25, 0.75 ];

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
%% save data
save d_BEM d_BEM
save gamma_rad_BEM rad
save gamma_tot_BEM tot
save nearfield ee
save nearfieldx ex
save nearfieldy ey
save nearfieldz ez
save Qx Qx
save Qy Qz
save Qz Qz
save Qavg Q_avg
% save FEx FEx
% save FEy FEy
% save FEz FEz
save FE FE
end 
