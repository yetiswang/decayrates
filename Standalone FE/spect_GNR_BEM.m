%% BEM retarded calculation of nanorod scattering spectrum
% The particle is set by default to lie along the x axis, and plane wave exctation is
% set by default polarized along the x axis. The long axis has to be the
% height of the nanorod.
%       Input :   height and diamter of the nanorods
%                 enei: wavelength for evaluation, generated with linspace
%                 ztab: can be any number, and will be set to 0 in
%                 function.
%       Output:   sca: scattering spectrum with the same length as enei,
%                 fit: fitted lorentz curve, same length as enei,
%                 Lorentz: four Lorentz parameters
function [ sca, abs, fit, Lorentz, p] = spect_GNR_BEM( epstab, height, diameter, enei, ztab )

if nargin == 5
    ztab = 0 ;
else
    ztab = 'NaN';
end
% %% initialization
% eps_m = 1.33^2 ;
% %  table of dielectric functions
% epstab = { epsconst( eps_m ), epstable( 'gold.dat' ), epstable( 'bk7.dat' ) };

%  nanorod meshes
mesh = [ 24, 24, 26];
if ztab == 'NaN'
    %  options for BEM simulations
    op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
    fprintf('Simulation started... Substrate is not included. \n');
    
else
    
    %  default options for layer structure
    op = layerstructure.options;
    %  set up layer structure
    layer = layerstructure( epstab, [ 1, 3 ], ztab, op );
    %  options for BEM simulations
    op = bemoptions( 'sim', 'ret', 'interp', 'curv' , 'layer', layer );
    fprintf('Simulation started... Substrate is included. \n');
    
end

ndir = sprintf('%d nm %d nm',height,diameter);

% particle gold nanorod
p = trirod( diameter, height, mesh );
%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);
% if ztab is set, elevate the particle to 2 nm above the ztab
if ztab == 'NaN'
else
    p = shift( p, [ 0, 0, - min( p.pos( :, 3 ) ) + 2 + ztab ] );
end

% %  plot particle boundary
% figure
% plot( p, 'EdgeColor', 'b');
% title (ndir)
% axis on;
% axis tight;
% grid on;

%  set up COMPARTICLE object
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  TM mode, excitation from above
dir = [ 0, 0, -1];
pol = [ 1, 0, 0 ];

if ztab == 'NaN'
    
else
    %%  tabulated Green functions
    %  For the retarded simulation we first have to set up a table for the
    %  calculation of the reflected Green function.  This part is usually slow
    %  and we thus compute GREENTAB only if it has not been computed before.
    if ~exist( 'greentab', 'var' ) || ~ greentab.ismember( layer, enei, p )
        %  automatic grid for tabulation
        tab = tabspace( layer, p );
        %  Green function table
        greentab = compgreentablayer( layer, tab );
        %  precompute Green function table
        %    for a more accurate simulation of the layer the number of
        %    wavelenghts should be increased
        greentab = set( greentab, linspace( 400, 1050, 5 ), op );
    end
    op.greentab = greentab;
end
%%  BEM solver
%  initialize BEM solvers
bem = bemsolver( p, op );

%  initialize plane wave excitations
exc = planewave( pol, dir, op );
%  scattering cross sections
sca = zeros( numel( enei ), 1 );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavelengths
for ien = 1 : length( enei )
    %  surface charges
    sig = bem \ exc( p, enei( ien ) );
    %  scattering cross sections
    sca( ien, : ) = exc.sca( sig );
    abs( ien, : ) = exc.abs( sig );
    multiWaitbar( 'BEM solver', ien / numel( enei ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%% LORENTZIAN FIT
[ fit, Lorentz, ~ ] = LorentzFitYW( sca, enei, ndir, 'sca') ;
%[ fit_abs, Lorentz_abs ] = LorentzFitYW( abs, enei, ndir, 'abs', 'Yes' ) ;