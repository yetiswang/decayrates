%% BEM calculation of spheres

clear all
close all
clc

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  table of dielectric functions
epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };
%  diameter of sphere
diameter = 140;
Q_int = 0.01 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 41 ),  ...
    pi * linspace( 0, 1, 41 ) .^ 2, diameter );

%  choose distance dependent calculations or wavelength dependent
%  calculations
distance = 1;
wavelength = 0;
%%
% Get current folder
dir = pwd;
%%
if distance == 1
    % make a new directory
    ndir = ['D',num2str(diameter),'_GNS_distance_dependent_decayrates'];
    mkdir(ndir)
    cd(ndir)
    %% Calculate and save the scattering spectrum for reference
    [ sca, fit, Lorentz ] = spect_GNS_BEM( epstab, diameter, linspace(300,800,50 ));
    %%
    %  initialize sphere
    p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
    
    %%  dipole oscillator
    enei = 626;
    
    %1D positions of dipole
    %choose the minimum distance by setting the values of x vector
    % x = reshape( linspace( 0.515, (50 + 0.5 * diameter )/diameter, 50 ) * diameter, [], 1 );
    z = reshape( logspace( log10(0.51 * diameter), log10((50 + 0.5 * diameter )), 100 ) , [], 1 );
    
    %  compoint
    pt = compoint( p, [ 0 * z, 0 * z, z ], 'mindist' , 1e-10  );
    %  dipole excitation , x ,y and x direction
    dip = dipole( pt, [ 0, 0, 1 ; 0 , 1 , 0; 1,0,0 ], op );
    figure
    bem_plot(p)
    hold on
    plot3( 0 * z, 0 * z, z ,'r*')
    hold off
    saveas(gcf,'dipole-sphere.fig')
    saveas(gcf,'dipole-sphere.png')
    %%  BEM simulation
    %  set up BEM solver
    bem = bemsolver( p, op );
    %  surface charge
    sig = bem \ dip( p, enei );
    %  total and radiative decay rate
    % (in units of the free-space decay rate rad0)
    [ tot, rad, rad0 ] = dip.decayrate( sig );
    %% Compare results with Mie theory
    % Distance from surface to molecule
    d_BEM = z - diameter/2;
    mie = miesolver( epstab{ 2 }, epstab{ 1 }, diameter, op );
    %  total and radiative decay rate
    [ tot_mie, rad_mie ] = mie.decayrate( enei, z );
    
    %% Calculation of quantum yields
    
    QYz = rad(:,1)./(tot(:,1) + (1 - Q_int)/Q_int );
    QYx = rad(:,3)./(tot(:,3) + (1 - Q_int)/Q_int );
    QYy = rad(:,2)./(tot(:,2) + (1 - Q_int)/Q_int );
    rad_average = (rad(:,1) + rad(:, 2) + rad(:, 3) )./3 ;
    tot_average = (tot(:,1) + tot(:, 2) + tot(:, 3) )./3 ;
    QYavg = rad_average./( tot_average + ( 1 - Q_int )/Q_int );
    
    QYperp_mie = rad_mie(:,2)./(tot_mie(:,2) + (1 - Q_int)/Q_int );
    QYpara_mie = rad_mie(:,1)./(tot_mie(:,1) + (1 - Q_int)/Q_int );
    rad_mie_avg = (rad_mie(:,2) + 2.* rad_mie(:,1))./3;
    tot_mie_avg = (tot_mie(:,2) + 2.* tot_mie(:,1))./3;
    QYavg_mie = rad_mie_avg./(tot_mie_avg + (1 - Q_int)/Q_int );
    
    %%  final plot
    figure
    semilogy( z ,tot, '-', z, tot_mie, ' r. '  );  hold on;
    semilogy( z, rad, 'o-', z, rad_mie, ' b. ' );
    
    xlim( [ min( z ), max( z ) ] );
    
    title( 'Total and radiaitve decay rate for dipole oriented along x, y and x' );
    legend( 'tot_{\perp} BEM','tot_{//} BEM','tot_{//} BEM', 'tot_{\perp} Mie','tot_{//} Mie','rad_{\perp} BEM','rad_{//} BEM','rad_{//} BEM','rad_{\perp} Mie','rad_{//} Mie'  )
    xlabel( 'Position (nm)' );
    ylabel( 'Decay rate' );
    saveas(gcf, 'decayrates_d.fig')
    saveas(gcf, 'decayrates_d.png')
    figure
    plot(d_BEM, QYx, d_BEM, QYy, d_BEM, QYz, d_BEM, QYavg, d_BEM, QYavg_mie)
    leg = legend('$\phi_{\perp} BEM$','$\phi_{//} BEM$','$\phi_{//} BEM$','$\overline{QY} BEM$','$\overline{QY} Mie$');
    set(leg, 'Interpreter','latex');
    set(leg, 'FontSize',14);
    xlabel('Distance to surface (nm)')
    ylabel(['Quantum yield (intrinsic QY = ', num2str(Q_int), ' )'])
    title(['Wavelength of dipole ', num2str(enei),' nm'])
    xlim( [ 0 , max(d_BEM) ]  )
    
    hline(Q_int, 'r:','Intinsic QY')
    saveas(gcf, 'QY_enhancement_d.fig')
    saveas(gcf, 'QY_enhancement_d.png')
    
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


cd(dir)
% make a new directory
ndir = ['D',num2str(diameter),'_GNS_distance_dependent_decayrates'];
mkdir(ndir)
cd(ndir)
if wavelength == 1
    
    %  initialize sphere
    p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
    
    %%  dipole oscillator
    enei = linspace( 400, 900, 50 );
    %  compoint
    d = 3 ;
    pt = compoint( p, [ 0, 0, 0.5 * diameter + d ] );
    %  dipole excitation
    dip = dipole( pt, [ 0, 0, 1 ; 0 , 1 , 0; 1, 0, 0 ], op );
    %  initialize total and radiative scattering rate
    [ tot_wl, rad_wl, rad0_wl] = deal( zeros( numel( enei ), 2 ) );
    
    figure
    bem_plot(p)
    hold on
    plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3)  ,'r*')
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
        [ tot_wl( ien, : ), rad_wl( ien, : ), rad0_wl( ien, : ) ] = dip.decayrate( sig );
        
        multiWaitbar( 'BEM solver', ien / numel( enei ) );
    end
    %  close waitbar
    multiWaitbar( 'CloseAll' )
    
    %% Calculate quantum yields
    QYz_wl = rad_wl(:,1)./(tot_wl(:,1) + (1 - Q_int)/Q_int );
    QYy_wl = rad_wl(:,2)./(tot_wl(:,2) + (1 - Q_int)/Q_int );
    QYx_wl = rad_wl(:,3)./(tot_wl(:,3) + (1 - Q_int)/Q_int );
    rad_wl_avg = (rad_wl(:, 1) + rad_wl(:, 2) + rad_wl(:, 3) )./3;
    tot_wl_avg = (tot_wl(:, 1) + tot_wl(:, 2) + tot_wl(:, 3) )./3;
    QY_avg_wl = rad_wl_avg./(tot_wl_avg +(1 - Q_int)/Q_int ) ;
    
    %%  final plot
    % quantum yield averaged over all solid angles
    figure
    plot( enei, tot_wl       );  hold on;
    plot( enei, rad_wl, 'o-' );
    xlabel( 'Wavelength (nm)' );
    ylabel( 'decay rates' );
    legend( 'tot_{\perp} BEM','tot_{//} BEM','tot_{//} BEM','rad_{\perp} BEM','rad_{//} BEM','rad_{//} BEM' )
    saveas(gcf, 'decayrates_wl.fig')
    saveas(gcf, 'decayrates_wl.png')
    figure
    plot(enei, QYz_wl , enei, QYy_wl, enei, QYs_wl, enei, QY_avg_wl )
    hline(Q_int, 'r--','Intrinsic QY')
    leg = legend('$\phi_{\perp} BEM$','$\phi_{//} BEM$','$\phi_{//} BEM$','$\overline{QY} BEM$');
    set(leg, 'Interpreter','latex');
    set(leg, 'FontSize',14);
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
