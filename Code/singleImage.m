%--------------GRABBING FILE DATA--------------------%
clear
clf 

filename = 'bunny_shift_5_5_5.xls';
bunny=(xlsread(filename)./5-1).*100+5; %Reads in the bunny data, shifted to first octant
sourcepool=nodePoolGenGrid(bunny,.5,1.5); %pool of nodes gridded
eval=bunny(1:300:34834,:); % Generate pool here so it doesn't need to be redone
colopool=bunny(5:300:34834,:); % Note that evals != colos; offset

%-------------SETTING UP PROBLEM GEOMETRY
n=5; % # of nodes
tic;
[heats,sources,colos]=bunnyGenError(n,sourcepool,eval,colopool,bunny); %Get error values
toc
scatter3(bunny(:,1),bunny(:,2),bunny(:,3),1,heats); %Generate bunny figure
hold on; %retain current plot when generating new plots
%scatter3(colos(:,1),colos(:,2),colos(:,3),20,[0,0,0], 'filled'); 

set(gcf,'color','w'); %background to clear

%---------------PLOTTING OPTIONS--------------%
% xLow=.9.*min(sourcepool(:,1));
% xHigh=1.1.*max(sourcepool(:,1));
% yLow=.9.*min(sourcepool(:,2));
% yHigh=1.1.*max(sourcepool(:,2)); %SHOWS SOURCES
% zLow=.9.*min(sourcepool(:,3));                    %UNCOMMENT 23-32 together
% zHigh=1.1.*max(sourcepool(:,3));
% xlim([xLow, xHigh]);
% ylim([yLow, yHigh]); 
% zlim([zLow, zHigh]);
% scatter3(sources(:,1),sources(:,2),sources(:,3),10,[0,0,0],'filled');

% str=cellstr(num2str([1:n]'));
% text(sources(:,1),sources(:,2),sources(:,3)-.01,str); %number the nodes

% colormap and bar
colormap(hsv); %OPTIONS: hsv, cool, parula, hot, jet, some lame ones
colorbar;
%caxis([0 .019]); %use if plotting error
%caxis([.075 .15]); %use if plotting approximation

