function [tot_div_time, tot_cur_time] = total_divcurl_time(mydir, dt, ntimesteps, ndt, nrods)

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

xs = unique(grid(:,1));
ys = unique(grid(:,2));
U = zeros(length(xs),length(ys));
V = zeros(length(xs),length(ys));

Mis = zeros(int8(ntimesteps/ndt), grid_index, 4);

for i = 0:ndt:ntimesteps-ndt
    
    sprintf(strcat('Interpolating frame ',num2str(i),'\n'))
    
    % Calculate rod positions at time i
    dat = importdata(dfile,'\t', (nrods+1)*i + 1, 'bufsize', 10000000);
    rods = dat.data;
    
    % Calculate rod positions at time i+1
    dat = importdata(dfile,'\t', (nrods+1)*(i+1) + 1, 'bufsize', 10000000);
    rodsdt = dat.data;
   
    dr = [rods(:,2) , rods(:,1) ,  rodsdt(:,2) , rodsdt(:,1) ];
    
    % run the vector interpolation script
    m = vectorFieldSparseInterpPatrick( dr, grid, sqrt(2), 1, []);
    
    x=0;

    for i=1:length(xs)
        for j = 1:length(ys)
            x = x + 1; %(i-1)*length(ys) + j;
            U(i,j) = m(x,4)-m(x,2);
            V(i,j) = m(x,3)-m(x,1);
        end
    end

    div = divergence(xs, ys, U, V);
    cur = curl(xs, ys, U, V);

    t=int8(i/ndt + 1);
    tot_div_time(t,:) = [(t-1)*dt, sum(sum(div(isfinite(div))))];
    tot_cur_time(t,:) = [(t-1)*dt, sum(sum(cur(isfinite(cur))))];
    
end
