function Mi=vectorFieldSparseInterpPatrick(M,Pg,threshold,d0,polygon)
% vectorFieldSparseInterp interpolates a vector field on a user-specified grid (using a sparse correlation matrix)
%
% For a vector field v:  
%    Interpolation with a correlation matrix K : <v> = K * v      (convolution)
%    where:  K = sG*exp(-(dx^2+dy^2)/d0^2), sG = weigth vector
%
% SYNOPSIS   Mi=vectorFieldSparseInterp(M,Pg,threshold,d0,polygon)
%
% INPUT      M         : vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                        (where (y0,x0) is the base and (y,x) is the tip of
%                        the vector).
%            Pg        : regular grid points, stored in a (mx2)-matrix of the form [yg xg]m.
%            threshold : radius to calculate the sparse distance matrix (see help
%                        createSparseDistanceMatrix)
%            d0        : parameter for the weight function G=exp(-D.^2/(1+d0^2)),
%                        where D is the distance matrix between all grid
%                        points and all vector (base) positions.
%                        d0 must be a scalar.
%            polygon   : (optional - pass polygon=[] to disable). The interpolated vector
%                        can be cropped to remove vectors outside a given region of interest.
%                        To create the polygon use the functions ROIPOLY or
%                        GETLINE. These functions return the polygon vertices
%                        stored in two vectors y and x. Set polygon=[x y] to
%                        use with vectorFieldInterp.
%
% OUTPUT     Mi        : interpolated vector field.
%
% DEPENDENCES            vectorFieldSparseInterp uses { createDistanceMatrix (C-MEX function) }
%                        vectorFieldSparseInterp is used by { }
%
% Aaron Ponti, 02/11/2004
% Edited Patrick Oakes 02/07/2011 - Made it work without the MEX file.  Probably not smart to try it on a
% really large matrix, but for an image it's fine. 

% Check d0
if size(d0)~=[1 1]
    error('The input parameter d0 must be a scalar.');
end

% Vector base positions
Pi=M(:,1:2);

% Vectors
V=[M(:,3)-M(:,1) M(:,4)-M(:,2)];

% Calculate distances
% D=createSparseDistanceMatrix(Pg,Pi,threshold); - this is the original MEX file line
newD = [];
for i = 1:length(Pg)
    pt = Pg(i,:);
    pt = repmat(pt,length(Pi),1);
    dist = sqrt((pt(:,1)-Pi(:,1)).^2+(pt(:,2)-Pi(:,2)).^2);
    keep = find(dist<threshold);
    if size(keep,1)>0
    compare = [ones(length(keep),1)*i keep dist(keep)];
    newD = vertcat(newD,compare);
    end
    clear compare keep dist
end
    
D = sparse(newD(:,1), newD(:,2),newD(:,3),length(Pg),length(Pi));
clear newD pt compare

% Correlation matrix (d0 may be a scalar or a matrix)
G=D; % Copy the sparse matrix
G(find(D))=exp(-D(find(D)).^2./d0.^2); clear D;

% Interpolate
Vi=[G*V(:,1) G*V(:,2)];

% Normalize
sG=sum(G,2);

%When the distances from a grid point to all the data points are >
%threshold, the corresponding 'sG' is zero. In this case, we set it to be
%'NaN' to indicate that this grid point is not interpolated.
sG(find(sG==0))=NaN; % It also prevents division by zero
Vi=[Vi(:,1)./sG Vi(:,2)./sG];

% Mi is the interpolated M
%Mi=[Pg Pg+Vi];
Mi=[Pg Vi];

% Set all vectors outside the passed polygon to 0
if ~isempty(polygon)
    index=inpolygon(Mi(:,1),Mi(:,2),polygon(:,2),polygon(:,1));
    Mi(find(~index),3)=Mi(find(~index),1);
    Mi(find(~index),4)=Mi(find(~index),2);
end
