function [grid, Mis]=getAllVecFields(mydir, ntimesteps, ndt, nrods)
basedir='/home/simonfreedman/Code/cytomod/shila/semiflexible/out/network/';
dfile  =strcat(basedir,mydir,'/txt_stack/rods.txt');
% make a grid
xrange = 50; %um
yrange = 50; %um
dt = 0.001;
grid = zeros(xrange*yrange, 2);
grid_index = 1;
for x = (-xrange/2):1:(xrange/2)
    for y = (-yrange/2):1:(yrange/2)
        grid( grid_index , : ) = [y,x];
        grid_index = grid_index + 1;
    end
end

Mis = zeros(int8(ntimesteps/ndt), grid_index - 1, 4);
% Import the file
dat = importdata(dfile,'\t', 1);
mat = dat.data;
angs = atan2(mat(:,3),mat(:,4));
cosangs = cos(angs);
sinangs = sin(angs);
rods = [mat(:,2) - sinangs, mat(:,1) - cosangs, mat(:,2) + sinangs, mat(:,1) + cosangs];

for i = 2:ntimesteps
    
    sprintf(strcat('Interpolating frame ',num2str(i),'\n'))
    
    % Calculate new rod positions
    dat = importdata(dfile,'\t', (nrods+1)*(i-1) + 1);
    mat = dat.data;
    angs = atan2(mat(:,3),mat(:,4));
    cosangs = cos(angs);
    sinangs = sin(angs);
    rodsdt = [mat(:,2) - sinangs, mat(:,1) - cosangs, mat(:,2) + sinangs, mat(:,1) + cosangs];
    
    % run the vector interpolation script
    m = vectorFieldSparseInterpPatrick(rodsdt - rods, grid, sqrt(2)/2, 1,[]);
    if mod(i, ndt) == 0
        quiver( grid(:,2) ,grid(:,1) , m(:,4)-m(:,2) , m(:,3)-m(:,1) );
        title(strcat('t = ', str2num(i*dt)));
        xlabel('\mu m')
        ylabel('\mu m')
        frames((i-1)/ndt + 1) = getframe(gcf);
    end

    Mis(int8((i-1)/ndt),:,:)=m;
    rods = rodsdt;

end

movie2avi(frames, strcat(basedir, mydir, '/data/velocity_vectors.avi'));
