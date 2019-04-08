function main_FE_1D
% 
%height = [  35 40 45 50 55 58 60 62 65 70 73 75 80 83 90];
%diameter = [ 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30];

%height = [ 70 75 80 85 90 100 105 108 115 120 135 140 150 155 160  ];
%diameter = [ 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 ];

height = [ 35 65 87 105 120 128 ];
diameter = [10 20 30 40 50 60];
%height = [  90 100 110 115 118 123 130 135 140 145 150 155];
%diameter = [ 45 45 45 45 45 45 45 45 45 45 45 45];

%height = [52 55 55 43 58 67 58 67 58 51 61 54 45 34 60 69 66 47 53 67 58 66 56 64 74 64];
%diameter = [28 30 30 29 25 31 30 31 24 29 28 26 21 28 31 30 24 25 31 32 31 31 28 27 31 31];
%height = [66 56 64 74 64];
%diameter = [31 28 27 31 31];

%height = 60 ;
%diameter = 30 ;

metal = {'AgJC'};
parent = pwd; 
for i = 1 : length(metal)
    for j = 1 : length(height)
        if isequal(height(j),diameter(j))
            FE_GNS_1D(diameter(j),metal{i});
            close all
        else
            [PCR{j},FE{j},Lorentz{j},ee{j}] = FE_GNR_1D(height(j),diameter(j),metal{i},637, 670);
            %FE{j} = FE_GNR_1D(height(j),diameter(j),metal{i});
            close all
        end
        cd(parent)
        MaximumPCRenhancement(j,:) = PCR{j}.maximumPCRenhancement;
        Isat(j) = PCR{j}.I_satenh(10)./ee{j}(10);
        WL(j) = 1248./Lorentz{j}(3);
    end
    cd(parent)
end

%% Plotting
figure
for i = 1 : 6
    %scatter(WL,MaximumPCRenhancement(:,i))
    %hold on
    plot(WL,MaximumPCRenhancement(:,i),'o-')
    hold on
    xlabel('SPR wavelength (nm)')
    ylabel('Maximal photon rate enhancement')
end
legend({'\phi_0 = 0.01','\phi_0 = 0.02','\phi_0 = 0.05','\phi_0 = 0.10','\phi_0 = 0.40','\phi_0 = 1.00'})
hold off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(ax,'YScale','log')
vline(637,'-','Exc.')
vline(670,'-','Dip.')
saveas( gcf,'Photon_rate enhancement.fig' )
saveas( gcf,'Photon_rate enhancement.png' )
%%
figure
%scatter(WL,MaximumPCRenhancement(:,i))
%hold on
semilogy(WL,Isat,'o-')
hold on
xlabel('SPR wavelength (nm)')
ylabel('Saturation intensity (W/m^2)')

hold off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(ax,'YScale','log')
hline(PCR{1,1}.I_sat)
saveas( gcf,'Isat enhancement.fig' )
saveas( gcf,'Isat enhancement.png' )
%% 
figure

plot(diameter,MaximumPCRenhancement(1:6,1),'o-')
hold on
xlabel('Diameter (nm)')
ylabel('Maximal photon rate enhancement')
title('QY = 0.01, SPR = 650 nm')
%legend({'\phi_0 = 0.01','\phi_0 = 0.02','\phi_0 = 0.05','\phi_0 = 0.10','\phi_0 = 0.40','\phi_0 = 1.00'})
hold off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(ax,'YScale','log')
saveas( gcf,'Photon_rate enhancement_volume_depend.fig' )
saveas( gcf,'Photon_rate enhancement_volume_depend.png' )
%% Saving results
save FE
save Lorentz
save I
save PCR
save WL
save ee
save MaixmumPCRenhancement
end