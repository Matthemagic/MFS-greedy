function [errorList, sourcef, colof]=bunnyGenApprox(n,sourcepool,eval,colopool,bunny)
%Performs greedy optimization algorithm for source and
%colocation points pairings

%BC=@(x,y,z)1/sqrt(x.^2+y.^2+z.^2)+5.*x-2.*y+3.*z.^2+15*sin(x); %Define the true boundary condition
%BC=@(x,y,z)x.^2-y.^2-z.^2;
BC=@(x,y,z)1./sqrt(x.^2+y.^2+z.^2);

evals=[eval(:,1) eval(:,2) eval(:,3) BC(eval(:,1),eval(:,2),eval(:,3))];

if n>size(colopool,1) || n>size(sourcepool,1)  %Sanity check
    n = min(size(colopool,1),size(sourcepool,1));
end

sourcef=zeros(n,3); %Preallocated for speed originally
colof=zeros(n,3); %DEALLOCATE THIS FOR ERRORFN, REWRITE L8R

bestc=1; %These are pointers for the row which contains the best node/point
bestc1=1;
bests=1;

minerrorC=0; %Stores minimum error for a given marriage
minerrorS=0;

error=0;

errorList=zeros(n,1);

for i=1:n
    for j=1:size(sourcepool,1)
        sourcef(i,:)=sourcepool(j,:);
        for k=1:size(colopool,1)
            colof(i,:)=colopool(k,:);
            error=errorFn(sourcef(1:i,:),colof(1:i,:),evals);
            if k==1
                minerrorC=error;
                bestc1=1;
            elseif error < minerrorC
                minerrorC=error;
                bestc1=k;
            end
            %determine best colocation point for n source points
        end
        if j==1
            minerrorS=minerrorC;
            bests=1;
        elseif minerrorS > minerrorC
            bests=j;
            bestc=bestc1;
            minerrorS=minerrorC;
        
        end
        %determine best n-node approximation
    end
    %Adds the best source and colocation to their respective lists
    %Then deletes them from the pool of potential candidates
    sourcef(i,:)=sourcepool(bests,:);
    sourcepool(bests,:)=[];
    colof(i,:)=colopool(bestc,:);
    colopool(bestc,:)=[];
    errorList(i)=minerrorS;
end

%Run the error function one more time to display the coefficients

BCs=BC(colof(:,1),colof(:,2),colof(:,3)); %

A=zeros(size(colof,1));
for j=1:size(colof,1)
    A(j,:)=1./sqrt((sourcef(:,1)-colof(j,1)).^2+(sourcef(:,2)-colof(j,2)).^2+(sourcef(:,3)-colof(j,3)).^2);
end

coeffs=linsolve(A,BCs);

F=@(x,y,z)0;
for j = 1:size(sourcef,1)
    l=@(x,y,z)coeffs(j)./sqrt((sourcef(j,1)-x).^2 + (sourcef(j,2)-y).^2 + (sourcef(j,3)-z).^2);
    F=@(x,y,z)(F(x,y,z)+l(x,y,z));
end

%Return approximation values
heats=F(bunny(:,1),bunny(:,2),bunny(:,3));

end