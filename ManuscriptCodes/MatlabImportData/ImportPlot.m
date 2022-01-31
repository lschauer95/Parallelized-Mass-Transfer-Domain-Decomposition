    
[figSettings, color, rgbcmy, mList, lineList, dashList] = setFigureDefaults();

% simulation parameter input
D=1;
dt=0.1;
Np = 10000000;
time = 10;

% number of time steps the user writes to file -- used to shape files
nWrites = 1;

spacemassplot = true;
import = false;

if(spacemassplot)
    
    
    if(import)
        locsXfileopen = fopen('locsX.dat');
        locsXfile = fscanf(locsXfileopen,'%f');
        locsXfile(isnan(locsXfile)) = [];
        locsXfileshaped = reshape(locsXfile,Np,1,nWrites);
    %     X = locsXfileshaped;

        locsYfileopen = fopen('locsY.dat');
        locsYfile = fscanf(locsYfileopen, '%f');
        locsYfile(isnan(locsYfile)) = [];
        locsYfileshaped = reshape(locsYfile,Np,1,nWrites);
    %     Y = locsYfileshaped;

        massfileopen = fopen('mass.dat');
        massfileraw = fscanf(massfileopen,'%f');
        massfileraw(isnan(massfileraw)) = [];
    %     massfileshaped = reshape(massfileraw,SqrtNp,SqrtNp,51);
        massvector = reshape(massfileraw,Np,1,nWrites);
    end

    analytical  = zeros(Np,1);
    
    % You can start a loop here over i with the final entry of x, y, z, and
    % mass being "i" to create a gif. This is currently hard coded to output
    % just a single time step. 
    
    x = locsXfileshaped(:,1,1);
    y = locsYfileshaped(:,1,1);
%     z = locsZfileshaped(:,1,4);
%         massvector = reshape(massfileshaped(:,:,i),[1,Np]);
    mass = massvector(:,1,1);
    
    xmidpt = 500;
    analytical = 0.5 * erfc(-((x - xmidpt)/ sqrt(4.0 * D * time)));
    
    subplot(2,2,[1 2])
    scatter(x, y, [], mass)
    hold on        
%         colormap(hot)
    c = colorbar();
    c.Location = 'northoutside';
    c.Label.String = 'Mass';
%     set(get(colorbar,'Title'),'String','Mass');
    set(gca,'xlim', [0, 1000], 'ylim', [0, 1000]);
    xlabel('X')
    ylabel('Y')
%     ylabel(colorbar,'Mass','Rotation',270);
%     hColourbar.Label.Position(1) = 5;
    hold off

    subplot(2,2,[3 4])
    scatter(x, mass)
    hold on
    scatter(x, analytical)
    xlabel('X')
    ylabel('Mass')
    legend('Observed Mass Values','Analytical Solution','Location','northwest')
    hold off

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     saveas(gcf, ['FinalTimeSolution'],'epsc')

end




