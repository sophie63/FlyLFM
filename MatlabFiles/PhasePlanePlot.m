[A,B,C]=pca(TSo);
folder = uigetdir();
figure
hold on
for i=1:20651
plot(B(i,1),B(i,2),'o','MarkerSize',1,'color',[0,0,Walk(i)/max(Walk)])

    if ~mod(i,5)
%M(i)=getframe(gcf);
    saveas(gcf,strcat(folder,'/100133_',num2str(i/5),'.tif'));
end
end