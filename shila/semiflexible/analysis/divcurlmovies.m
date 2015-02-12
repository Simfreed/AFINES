function divCurlMovies(mydir, ntimesteps, ndt, nrods, tfinal)

basedir='/home/simonfreedman/Code/cytomod/shila/semiflexible/out/network/';

[g,vs] = getVecFields(mydir, ntimesteps, ndt, nrods);
s = size(vs);

for i = 1:s(1,1)

    [x, y, d, c] = divcurl(g, vs(i,:,:));
    imagesc(x,y,d);
    title(strcat('t = ',num2str(tfinal*i/s(1,1))));
    colorbar;
    xlabel('\mu m');
    ylabel('\mu m');
    dframes(i) = getframe(gcf);

    imagesc(x,y,c);
    title(strcat('t = ',num2str(tfinal*i/s(1,1))));
    colorbar;
    xlabel('\mu m');
    ylabel('\mu m');
    cframes(i) = getframe(gcf);

end

movie2avi(cframes, strcat(basedir, mydir, '/data/curl.avi'));
movie2avi(dframes, strcat(basedir, mydir, '/data/div.avi'));

