function [grid, Mis]=getVecFields(mydir, ntimesteps, ndt, nrods, toPlot)
basedir='/home/simonfreedman/Code/cytomod/shila/semiflexible/out/network/';
dfile  =strcat(basedir,mydir,'/txt_stack/rods.txt');
% make a grid
xrange = 50; %um
yrange = 50; %um
dt = 0.001;
grid = zeros(xrange*yrange, 2);
grid_index = 0;
for x = (-xrange/2):1:(xrange/2)
    for y = (-yrange/2):1:(yrange/2)
        grid_index = grid_index + 1;
        grid( grid_index , : ) = [y,x];
    end
end

Mis = zeros(int8(ntimesteps/ndt), grid_index, 4);

for i = 0:ndt:ntimesteps-ndt
    
    sprintf(strcat('Interpolating frame ',num2str(i),'\n'))
    
    % Calculate rod positions at time i
    dat = importdata(dfile,'\t', (nrods+1)*i + 1);
    rods = dat.data;
    
    % Calculate rod positions at time i+1
    dat = importdata(dfile,'\t', (nrods+1)*(i+1) + 1);
    rodsdt = dat.data;
   
    dr = [rods(:,2) , rods(:,1) ,  rodsdt(:,2) , rodsdt(:,1) ];
    % run the vector interpolation script
    m = vectorFieldSparseInterpPatrick( dr, grid, sqrt(2), 1, []);
    
    if toPlot
        % Plot the vector field
        quiver( grid(:,2) ,grid(:,1) , m(:,4)-m(:,2) , m(:,3)-m(:,1) );
        title(strcat('t = ', num2str((i+1)*dt)));
        xlabel('\mu m')
        ylabel('\mu m')
        frames(int8(i/ndt + 1)) = getframe(gcf);

    end
    
    % Save the velocity vectors
    Mis(int8(i/ndt + 1),:,:)=m;
    
end

if toPlot
    movie2avi(frames, strcat(basedir, mydir, '/data/velocity_vectors.avi'));
end
