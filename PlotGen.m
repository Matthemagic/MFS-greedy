clf
a=semilogy(1:6,errors,'linewidth',2)
hold on

errorT=0.003278888.*ones(6,1)
b=semilogy(1:6,errorT,'linewidth',2)
set(gcf,'color','w'); %background to clear
xlabel('# of Nodes')
ylabel('Maximum Error')
legend([a b], 'Optimized Algorithm','Traditional Algorithm w/ n=6')