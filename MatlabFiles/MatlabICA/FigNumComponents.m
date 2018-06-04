figure1 = figure;
xtic=0.5+1:3:36;
% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'OL','VLNP','VMNP','AL','MB','LH','SNP','CX','LX','INP','PENP','GNG'},...
    'XTick',xtic,...
    'LineStyleOrderIndex',1);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 13]);
box(axes1,'on');
hold(axes1,'on');


plot((1:3:39)',GCaMP6,'o','MarkerFaceColor','b','MarkerSize',3,'color','b','Parent',axes1)
%plot((1:3:39)',GCaMP640x,'o','MarkerFaceColor','g','MarkerSize',3,'color','g')
plot((2:3:39)',Arclight,'o','MarkerFaceColor','r','MarkerSize',3,'color','r')
plot((1:3:39)',GCaMP6med,'s','MarkerFaceColor','b','color','b','MarkerSize',8)
plot((2:3:39)',Arclightmed,'s','MarkerFaceColor','r','color','r','MarkerSize',8)