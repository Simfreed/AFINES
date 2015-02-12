function [tot_div_time, tot_cur_time] = total_divergence_time(mdir, ntimesteps,ndt, nrods)

[g, velocities]=getVecFields(mdir, ntimesteps, ndt, nrods, 0);
xs = unique(g(:,1));
ys = unique(g(:,2));
U = zeros(length(xs),length(ys));
V = zeros(length(xs),length(ys));
tot_div_time = zeros(size(velocities,1), 2);

dt=0.2;
for t=1:size(velocities,1)
    x=0;

    for i=1:length(xs)
        for j = 1:length(ys)
            x = x + 1; %(i-1)*length(ys) + j;
            U(i,j) = velocities(t,x,4)-velocities(t,x,2);
            V(i,j) = velocities(t,x,3)-velocities(t,x,1);
        end
    end

    div = divergence(xs, ys, U, V);
    cur = curl(xs, ys, U, V);
    %size(div)
    tot_div_time(t,:) = [(t-1)*dt, sum(sum(div(isfinite(div))))];
    tot_cur_time(t,:) = [(t-1)*dt, sum(sum(cur(isfinite(cur))))];

end
