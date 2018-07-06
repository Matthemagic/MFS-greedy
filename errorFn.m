function error=errorFn(nodes,colPts,evals)
%Takes the running node and collocation points list

%Generates the approximation function and returns max error
U=@(x,y,z)1/sqrt(x.^2+y.^2+z.^2);
%U=@(x,y,z)x.^2-y.^2-z.^2;

BCs=zeros(size(colPts,1),1);
for j=1:size(colPts,1)
    BCs(j)=U(colPts(j,1),colPts(j,2),colPts(j,3));
end

%BCs=U(colPts(:,1),colPts(:,2),colPts(:,3)); %THIS VECTORIZATION NO WORKS

A=zeros(size(colPts,1));
for j=1:size(colPts,1)
    row=zeros(1,size(nodes,1));
    for h=1:size(nodes,1)
        row(h)=1/sqrt((nodes(h,1)-colPts(j,1)).^2+(nodes(h,2)-colPts(j,2)).^2+(nodes(h,3)-colPts(j,3)).^2);
    end
    A(j,:)=row;
end

%for j=1:size(colPts,1)
%   A(j,:)=1/sqrt((nodes(:,1)-colPts(j,1)).^2+(nodes(:,2)-colPts(j,2)).^2+(nodes(:,3)-colPts(j,3)).^2);
%end

coeffs=linsolve(A,BCs);

error=0;
for i = 1:size(evals,1)
    F=0;
    for j = 1:size(nodes,1)
        l=coeffs(j)./sqrt((nodes(j,1)-evals(i,1)).^2 + (nodes(j,2)-evals(i,2)).^2 + (nodes(j,3)-evals(i,3)).^2);
        F=F+l;
    end
    if abs(evals(i,4)-F) > error
       error=abs(evals(i,4)-F);
    end
end
end