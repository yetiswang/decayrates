%% BEM calculation of nanorods

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
height = [ 120  ] ;
diameter = [ 38  ];

Q_int = 0.9 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trirod ( diameter, height, mesh, 'triangles' );
%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);

    
%  choose distance dependent calculations or wavelength dependent
%  calculations
distance = 1;
wavelength = 0;

% Get current folder
dir = pwd;

if distance == 1
    % make a new directory 
    ndir = ['H',num2str(height),'D',num2str(diameter),'_GNR_distance_dependent_decayrates'];
    mkdir(ndir)
    cd( ndir )
    
    %%  dipole oscillator
    enei = 565;
    %1D positions of dipole
    %choose the minimum distance by setting the values of x vector
%     x = reshape( logspace( log10(0.51 * height), log10( 50 + 0.5 * height ), 100 ), [], 1 );
%     % compoint
%     pt = compoint( p, [ x, 0 .* x, 0 .* x ], 'mindist' , 1e-10  );
    z = reshape( logspace( log10(0.51 * diameter), log10( 50 + 0.5 * diameter ), 100 ), [], 1 ) ;
    pt = compoint( p, [ 0.*z, 0 .*z, z ], 'mindist' , 1e-10  );

    % dipole excitation , x and z direction
    dip = dipole( pt, [ 1, 0, 0 ; 0 , 1 , 0; 0, 0 ,1 ], op );
    
    figure
    bem_plot(p)
    hold on
    plot3(pt.pos(:,1), pt.pos(:,2),pt.pos(:,3) ,'r.')
    hold off
    
    %%  BEM simulation
    %  set up BEM solver
    bem = bemsolver( p, op );
    %  surface charge
    sig = bem \ dip( p, enei );
    %  total and radiative decay rate
    [ tot, rad, rad0 ] = dip.decayrate( sig );
    
    %% Calculate Quantum yields
    % Distance from surface to molecule
    %d_BEM = x - height/2;
    d_BEM = z - diameter/2;
    QYx = rad(:,1)./(tot(:,1) + (1 - Q_int)/Q_int );
    QYy = rad(:,2)./(tot(:,2) + (1 - Q_int)/Q_int );
    QYz = rad(:,3)./(tot(:,3) + (1 - Q_int)/Q_int );
    rad_average = (rad(:,1) + rad(:, 2) + rad(:, 3))./3 ;
    tot_average = (tot(:,1) + tot(:, 2) + tot(:, 3))./3;
    QYavg = rad_average./( tot_average + ( 1 - Q_int )/Q_int );
    
    %%  final plot
    figure
%     semilogy( x , tot, '-'  );  hold on;
%     semilogy( x, rad, 'o-' );
%     xlim( [ min( x ), max( x ) ] );
semilogy( z , tot, '-'  );  hold on;
semilogy( z, rad, 'o-' );
xlim( [ min( z ), max( z ) ] );

    
    title( 'Total and radiaitve decay rate' );
    legend( {'\gamma_{\perp}^{tot}/\gamma_0 BEM','\gamma_{||}^{tot}/\gamma_0 BEM', '\gamma_{||}^{tot}/\gamma_0 BEM','\gamma_{\perp}^{rad}/\gamma_0 BEM','\gamma_{||}^{rad}/\gamma_0 BEM','\gamma_{||}^{rad}/\gamma_0 BEM'} )
    xlabel( 'Position (nm)' );
    ylabel( 'Decay rate' );
    saveas(gcf, 'decayrates_wl.fig')
    saveas(gcf, 'decayrates_wl.png')
    
    figure
    plot(d_BEM, QYx, d_BEM, QYy, d_BEM, QYz, d_BEM , QYavg )
    legend({'$\phi_{\perp} BEM$', '$\phi_{||} BEM$','$\phi_{||} BEM$','$\overline{\phi} BEM$'},'interpreter','latex')
    xlabel('Distance to surface (nm)')
    ylabel(['Quantum yield (intrinsic QY = ', num2str(Q_int), ' )'])
    title(['Wavelength of dipole ', num2str(enei),' nm'])
    xlim( [ 0 , max(d_BEM) ]  )
    saveas(gcf, 'QY_enhancement_wl.fig')
    saveas(gcf, 'QY_enhancement_wl.png')
    %% save data
    save d_BEM d_BEM
    save gamma_rad_BEM rad
    save gamma_tot_BEM tot
    save QY_average_BEM QYavg
    save QY_x QYx
    save QY_y QYy
    save QY_z QYz
else
end


%% Wavelength dependent calculation
cd(dir)
% make a new directory
ndir = ['H',num2str(height),'D',num2str(diameter),'_GNR_wavelength_dependent_decayrates'];
mkdir(ndir)
cd(ndir)

if wavelength == 1
    %%  dipole oscillator
    enei = linspace( 400, 1000, 30 );
    %  compoint
    d = 7.2 ;
    %pt = compoint( p, [ (0.5 * height) + d, 0, 0 ] );
    pt = compoint( p, [ 0, 0, ( 0.5*diameter ) + d ] );
    %  dipole excitation
    dip = dipole( pt, [ 1, 0, 0; 0, 1, 0 ; 0, 0, 1 ], op );
    %  initialize total and radiative scattering rate
    [ tot_wl, rad_wl ] = deal( zeros( numel( enei ), 3 ) );
    
    figure
    bem_plot(p)
    hold on
    plot3(pt.pos(:,1), pt.pos(:,2),pt.pos(:,3) ,'r.')
    hold off
    %%  BEM simulation
    %  set up BEM solver
    bem = bemsolver( p, op );
    
    multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
    %  loop over wavelengths
    for ien = 1 : length( enei )
        %  surface charge
        sig = bem \ dip( p, enei( ien ) );
        %  total and radiative decay rate
        [ tot_wl( ien, : ), rad_wl( ien, : ) ] = dip.decayrate( sig );
        
        multiWaitbar( 'BEM solver', ien / numel( enei ) );
    end
    %  close waitbar
    multiWaitbar( 'CloseAll' );
    
    
    %% Calculate quantum yields 
    QYx_wl = rad_wl(:,1)./(tot_wl(:,1) + (1 - Q_int)/Q_int );
    QYy_wl = rad_wl(:,2)./(tot_wl(:,2) + (1 - Q_int)/Q_int );
    QYz_wl = rad_wl(:,3)./(tot_wl(:,3) + (1 - Q_int)/Q_int );
    rad_wl_avg = (rad_wl(:,1) + rad_wl(:,2) + rad_wl(:,3))/3;
    tot_wl_avg = (tot_wl(:,1) + tot_wl(:,2) + tot_wl(:,3))/3;
    QY_avg_wl = rad_wl_avg./(tot_wl_avg + (1 - Q_int)/Q_int ) ;

    %%  final plot
    % quantum yield averaged over all solid angles
    figure
    plot( enei, tot_wl       );  hold on;
    plot( enei, rad_wl, 'o-' );
    xlabel( 'Wavelength (nm)' );
    ylabel( 'decay rates' );
    legend( {'\gamma_{\perp}^{tot}/\gamma_0','\gamma_{||}^{tot}/\gamma_0', '\gamma_{||}^{tot}/\gamma_0','\gamma_{\perp}^{rad}/\gamma_0','\gamma_{||}^{rad}/\gamma_0','\gamma_{||}^{rad}/\gamma_0'} )
    saveas(gcf, 'decayrates_wl.fig')
    saveas(gcf, 'decayrates_wl.png')
    figure
    plot( enei, tot_wl_avg       );  hold on;
    plot( enei, rad_wl_avg, 'o-' );
    xlabel( 'Wavelength (nm)' );
    ylabel( 'decay rates' );
    legend( {'$\overline{\gamma}^{tot}/\gamma_0$','$\overline{\gamma}^{rad}/\gamma_0$'},'interpreter','latex' )
    saveas(gcf, 'decayrates_wl_avg.fig')
    saveas(gcf, 'decayrates_wl_avg.png')
    figure
    plot(enei, QYx_wl , enei, QYy_wl, enei, QYz_wl,enei, QY_avg_wl )
    hline(Q_int, 'r--','Intrinsic QY')
    legend({'$\phi_{\perp} $', '$\phi_{||} $','$\phi_{||} $','$\overline{\phi} $'},'interpreter','latex')
    xlabel('Wavelength (nm)')
    ylabel(['Quantum yield (intrinsic QY = ', num2str(Q_int), ' )'])
    xlim( [ min(enei) , max(enei) ]  )
    saveas(gcf, 'QY_enhancement_wl.fig')
    saveas(gcf, 'QY_enhancement_wl.png')
    %% save data
    save wavelength enei
    save gamma_rad_BEM rad_wl
    save gamma_tot_BEM tot_wl
    save QY_average_BEM QY_avg_wl
    save QYx_wl QYx_wl
    save QYy_wl QYy_wl
    save QYz_wl QYz_wl
else
end
