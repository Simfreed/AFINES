function [xs,ys,div,cur] = divcurl(grid, ms) %mydir, ntimesteps, ndt, nrods)
%[g, ms]=getVecFields(mydir, ntimesteps, ndt, nrods)
xs = unique(grid(:,1));
ys = unique(grid(:,2));
U = zeros(length(xs),length(ys));
V = zeros(length(xs),length(ys));
x=0;
for i=1:length(xs)
    for j = 1:length(ys)
        x = x + 1; %(i-1)*length(ys) + j;
        U(i,j) = ms(1,x,4)-ms(1,x,2);
        V(i,j) = ms(1,x,3)-ms(1,x,1);
    end
end
div = divergence(xs, ys, U, V);
cur = curl(xs,ys,U,V);
