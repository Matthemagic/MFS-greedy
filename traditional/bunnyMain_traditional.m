%Main file. Performs greedy optimization algorithm for source and
%colocation points pairings
%Dependencies: gensurfacesphere, gensourcesphere, errorFunctionsphere
warning off
format long
clear;
clf;
tic;
n=6; %Number of final colocation-source marriages
BC=@(x,y,z)1./sqrt(x.^2+y.^2+z.^2)+5; %Define the true boundary condition
%BC=@(x,y,z)x.^2-y.^2-z.^2;

filename = 'bunny_shift_5_5_5.xls';
bunny=(xlsread(filename)./5-1).*100+5;
eval=bunny(1:1500:34834,:);
%eval=nodePoolGenGrid(bunny,.5,1.6); %pool of 441 nodes gridded

colopool=bunny(5:6714:34834,:);
nodes=[mean(bunny(:,1)),mean(bunny(:,2)),(mean(bunny(:,3))+range(bunny(:,3)));
    mean(bunny(:,1)),mean(bunny(:,2)),(mean(bunny(:,3))-range(bunny(:,3)));
    mean(bunny(:,1)),mean(bunny(:,2))+range(bunny(:,3)),(mean(bunny(:,3)));
    mean(bunny(:,1)),mean(bunny(:,2))-range(bunny(:,3)),(mean(bunny(:,3)));
    mean(bunny(:,1))+range(bunny(:,3)),mean(bunny(:,2)),(mean(bunny(:,3)));
    mean(bunny(:,1))-range(bunny(:,3)),mean(bunny(:,2)),(mean(bunny(:,3)))];
evals=zeros(size(eval,1),4);
for i=1:size(eval,1)
    evals(i,:)=[eval(i,1) eval(i,2) eval(i,3) BC(eval(i,1),eval(i,2),eval(i,3))];
end

colof=colopool; %DEALLOCATE THIS FOR ERRORFN, REWRITE L8R

errorList=zeros(n,1);

BCs=zeros(size(colof,1),1);
for j=1:size(colof,1)
    BCs(j)=BC(colof(j,1),colof(j,2),colof(j,3));
end

A=zeros(size(colof,1));
for j=1:size(colof,1)
    row=zeros(1,n);
    for h=1:n
        row(h)=1/sqrt((nodes(h,1)-colof(j,1))^2+(nodes(h,2)-colof(j,2))^2+(nodes(h,3)-colof(j,3))^2);
    end
    A(j,:)=row;
end
coeffs=linsolve(A,BCs);
error=0;
for i = 1:size(evals,1)
    F=0;
    for j = 1:n
        l=coeffs(j)./sqrt((nodes(j,1)-evals(i,1)).^2 + (nodes(j,2)-evals(i,2)).^2 + (nodes(j,3)-evals(i,3)).^2);
        F=F+l;
    end
    if abs(evals(i,4)-F) > error
       error=abs(evals(i,4)-F);
    end
    errorList(i)=error;
end
F=@(x,y,z)0;
for j = 1:size(nodes,1)
    l=@(x,y,z)coeffs(j)./sqrt((nodes(j,1)-x).^2 + (nodes(j,2)-y).^2 + (nodes(j,3)-z).^2);
    F=@(x,y,z)(F(x,y,z)+l(x,y,z));
end

disp('Collocation Points:')
disp(colof)
disp('Corresponding Nodes:')
disp(nodes)
disp('Corresponding Coefficients:')
disp(coeffs)
disp('Error Convergence:')
disp(error)

heat=F(bunny(:,1),bunny(:,2),bunny(:,3));
%heat=abs(F(bunny(:,1),bunny(:,2),bunny(:,3))-BC(bunny(:,1),bunny(:,2),bunny(:,3)));

scatter3(bunny(:,1),bunny(:,2),bunny(:,3),1,heat(:))

c = flag(3);
colormap(c);

ax.GridColor = 'blue';

set(gca,'Color','k');

%caxis([0 .019]); %use if plotting error
%caxis([.075 .15]); %use if plotting approximation

hold on
%s=scatter3(colof(:,1),colof(:,2),colof(:,3),20,[0,0,0],'filled');

set(gcf,'color','k'); %background to clear
%l=legend([s],'Optimized Collocation Points','Location','northwest');

% str=cellstr(num2str([1:n]'));
% t=text(nodes(:,1)-.1,nodes(:,2)+.15,nodes(:,3)-.18,str); %number the nodes
% set(gcf,'color','w'); %background to clear

% l=legend([s],'Sources','Location','northwest');

% xLow=.9.*min(nodes(:,1));
% xHigh=1.1.*max(nodes(:,1));
% yLow=.9.*min(nodes(:,2));
% yHigh=1.1.*max(nodes(:,2)); %SHOWS SOURCES
% zLow=.9.*min(nodes(:,3));                    %UNCOMMENT 23-32 together
% zHigh=1.1.*max(nodes(:,3));
% xlim([xLow, xHigh]);
% ylim([yLow, yHigh]); 
% zlim([zLow, zHigh]);

toc