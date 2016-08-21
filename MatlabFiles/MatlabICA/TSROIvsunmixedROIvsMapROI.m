% Zscore maps
for i=1:Npc
    GM1vn(:,i)=GM(:,i)/(sqrt(var(GM(:,i))));
end

GMz=GM;
GMz(GM1vn<2.5)=0;

% Remove bad PCs
Id=[1,3,5];
Rnew=R-u(:,Id)*s(Id,Id)*v(:,Id)';
for j=1:Npc
parfor i=1:S1(4)
TSzmap(i,j)=mean(Rnew(i,:).*GMz(:,j)');
end
j
end

TSzmapo=TSzmap(:,Order(Npc:-1:1));
save(strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'TSzmap'),'TSzmapo')

figure
hold on
for i=1:150
    plot(TSo(:,i)/sqrt(var(TSo(:,i)))+i*10,'r')
    plot(TSzmapo(:,i)/sqrt(var(TSzmapo(:,i)))+i*10+5,'b')
end
% 
% %GMzbin=GMz;
% %GMzbin(GMz~=0)=1;
% 
% % % Averaging in z scored binarized map
% % for j=1:Npc
% % parfor i=1:S1(4)
% % Rm(i,:,j)=R(i,:).*GMzbin(:,j)';
% % end
% % j
% % end
% % 
% % TSroi=sum(Rm,2);
% 
% % Just sum in maps
% 
% for j=1:Npc
% parfor i=1:S1(4)
% Rmap(i,:,j)=R(i,:).*GMz(:,j)';
% end
% j
% end
% TSzmap=sum(Rmap,2);
% % 
% % % Masked and unmixed
% % for j=1:Npc    
% %     TSmum(:,j)=squeeze(Rm(:,:,j)*pinv(icasig(j,:)));
% % end
% 
% % Weighted average by zscored map and unmixed
% for j=1:Npc
% TSzum(:,j)=Rmap(:,:,j)*pinv(GMz(j,:));
% end
% 
% 
% n=30
% figure 
% plot(TS(:,n)/sqrt(var(TS(:,n))),'b')
% hold on
% 
% plot(TSroi(:,1,n)/sqrt(var(TSroi(:,1,n))),'r')
% plot(TSmum(:,n)/sqrt(var(TSmum(:,n))),'g')
% plot(TSzmap(:,1,n)/sqrt(var(TSroi(:,1,n))),'r')
% %plot(TSmum(:,n)/sqrt(var(TSmum(:,n))),'g')
% 
% plot(TSzum(:,n)/sqrt(var(TSzum(:,n))),'m')