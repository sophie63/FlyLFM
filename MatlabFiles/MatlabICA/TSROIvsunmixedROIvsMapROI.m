% Zscore maps
for i=1:Npc
    GM1vn(:,i)=GM(:,i)/(sqrt(var(GM(:,i))));
end

GMz=GM;
GMz(GM1vn<2.5)=0;
<<<<<<< HEAD
GMzbin=GMz;
GMzbin(GMz~=0)=1;

% % Averaging in z scored binarized map
% for j=1:Npc
% parfor i=1:S1(4)
% Rm(i,:,j)=R(i,:).*GMzbin(:,j)';
% end
% j
% end
% 
% TSroi=sum(Rm,2);

% Just sum in maps
||||||| merged common ancestors
GMzbin=GMz;
GMzbin(GMz~=0)=1;

for j=1:Npc
parfor i=1:S1(4)
Rm(i,:,j)=R(i,:).*GMzbin(:,j)';
end
j
end

TSroi=sum(Rm,2);

=======
%GMzbin=GMz;
%GMzbin(GMz~=0)=1;

% % Averaging in z scored binarized map
% for j=1:Npc
% parfor i=1:S1(4)
% Rm(i,:,j)=R(i,:).*GMzbin(:,j)';
% end
% j
% end
% 
% TSroi=sum(Rm,2);

% Just sum in maps
>>>>>>> 57759d29a73933bc92c6aaedf808d46b8fabe631
for j=1:Npc
parfor i=1:S1(4)
Rmap(i,:,j)=R(i,:).*GMz(:,j)';
end
j
end
TSzmap=sum(Rmap,2);
% 
% % Masked and unmixed
% for j=1:Npc    
%     TSmum(:,j)=squeeze(Rm(:,:,j)*pinv(icasig(j,:)));
% end

% Weighted average by zscored map and unmixed
for j=1:Npc
TSzum(:,j)=Rmap(:,:,j)*pinv(GMz(j,:));
end


n=30
figure 
plot(TS(:,n)/sqrt(var(TS(:,n))),'b')
hold on
<<<<<<< HEAD
%plot(TSroi(:,1,n)/sqrt(var(TSroi(:,1,n))),'r')
%plot(TSmum(:,n)/sqrt(var(TSmum(:,n))),'g')
||||||| merged common ancestors
plot(TSroi(:,1,n)/sqrt(var(TSroi(:,1,n))),'r')
plot(TSmum(:,n)/sqrt(var(TSmum(:,n))),'g')
=======
plot(TSzmap(:,1,n)/sqrt(var(TSroi(:,1,n))),'r')
%plot(TSmum(:,n)/sqrt(var(TSmum(:,n))),'g')
>>>>>>> 57759d29a73933bc92c6aaedf808d46b8fabe631
plot(TSzum(:,n)/sqrt(var(TSzum(:,n))),'m')