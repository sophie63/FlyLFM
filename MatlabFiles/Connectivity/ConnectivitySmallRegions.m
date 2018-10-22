%TS=TS';
% Load data
for i=1:75
TS1(:,i)=TS(:,i)/sqrt(var(TS(:,i)));
end

TS=TS1;
TS=TS';
Acov=TS*TS';
% Ainv=pinv(Acov);
DTS1=diff(TS');
for i=1:75
DTS(:,i)=DTS1(:,i)/sqrt(var(DTS1(:,i)));
end

DTS2=[DTS' DTS(end,:)'];

Apart=DTS2*TS';
% Afinal=Apart*Ainv;
Adiff=DTS2*DTS2';
imagesc(Acov)


W=Acov;
n  = size(W,1);             % number of nodes
M  = 1:n;                   % initial community affiliations
Q0 = -1; Q1 = 0;            % initialize modularity values
while Q1-Q0>1e-5;           % while modularity increases
Q0 = Q1;                % perform community detection
%[M, Q1] = community_louvain(W, [], M,'negative_sym');
[M, Q1] = community_louvain(W, [], M,'negative_asym');
end

[Or Amod]=reorder_mod(W,M);
f=figure
imagesc(Amod)

