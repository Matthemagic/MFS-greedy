function nodePool = nodePoolGenRand(surfface, m, dist)
%INPUTS: surface, a list of points that define the surface
%           -> must be nx3
%        dist, max distance the box should be from the surface
%        m, the number of points to be distributed in the box initially
%OUTPUT: nodePool, a list of candidate locations for node placement
%
%Randomly distributes points around surface
%Creates logical array to index points in/on surface
%Deletes points that are in/on the surface
%
tess = convhulln(surfface); %Generates the convex hull of the surface
[nt,~] = size(tess);

maxX=max(surfface(:,1))+dist;
maxY=max(surfface(:,2))+dist;
maxZ=max(surfface(:,3))+dist;

minX=min(surfface(:,1))-dist;
minY=min(surfface(:,2))-dist;
minZ=min(surfface(:,3))-dist;

rng('default'); %OPTIONAL for replicability

%Generate m random points on [0,1]x[0,1]x[0,1] then normalize to objective
randpool=rand(m,3);
nodePool=[randpool(:,1).*(maxX-minX)+minX randpool(:,2).*(maxY-minY)+minY randpool(:,3).*(maxZ-minZ)+minZ];

% use vectorized cross product since we're in 3-d
ab = surfface(tess(:,1),:) - surfface(tess(:,2),:);
ac = surfface(tess(:,1),:) - surfface(tess(:,3),:);
nrmls = cross(ab,ac,2);
degenflag = false(nt,1);

nrmllen = sqrt(sum(nrmls.^2,2));
nrmls = nrmls.*repmat(1./nrmllen,1,3);
% center point in the hull
center = mean(surfface,1);
% any point in the plane of each simplex in the convex hull
a = surfface(tess(~degenflag,1),:);
dp = sum((repmat(center,nt,1) - a).*nrmls,2);
k = dp<0;
nrmls(k,:) = -nrmls(k,:);

% We want to test if:  dot((x - a),N) >= 0
% If so for all faces of the hull, then x is inside
% the hull.
aN = sum(nrmls.*a,2);

in = false(m,1); %Initialize logical index

% if m is too large, we need to worry about the
% dot product grabbing huge chunks of memory.
memblock = 1e6;
blocks = max(1,floor(m/(memblock/nt)));
aNr = repmat(aN,1,length(1:blocks:m));
for i = 1:blocks
   j = i:blocks:m;
   if size(aNr,2) ~= length(j)
      aNr = repmat(aN,1,length(j));
   end
   in(j) = all((nrmls*nodePool(j,:)' - aNr) >= 0,1)';
end

nodePool(in,:)=[]; %Delet the nodes w/in the convex hull

end