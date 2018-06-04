clear

load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/973/973Registration/973TS.mat')
load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/973/973Registration/973TSnorm.mat')
for i=1:12
TS(:,i)=TS(:,i)/TSnorm(i);
end
TS=TS';
figure
plotTS
RegionsLightArclight(:,1)=mean(TS(:,2200:(2200+50)),2)-mean(TS(:,(2200-50):2200),2);
RegionsOdorArclight(:,1)=mean(TS(:,7459:(7459+50)),2)-mean(TS(:,(7459-50):7459),2);
figure
plot(RegionsLightArclight)
clear TS
load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/986/986Registration/986TS.mat')
load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/986/986Registration/986TSnorm.mat')
for i=1:12
TS(:,i)=TS(:,i)/TSnorm(i);
end
TS=TS';
figure
plotTS
RegionsLightArclight(:,2)=mean(TS(:,1000:(1000+50)),2)-mean(TS(:,(1000-50):1000),2);
RegionsOdorArclight(:,2)=mean(TS(:,7500:(7500+50)),2)-mean(TS(:,(7500-50):7500),2);
load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/100051/100051Registration/100051TS.mat')
load('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/ArclightCombo/100051/100051Registration/100051TSnorm.mat')
for i=1:12
TS(:,i)=TS(:,i)/TSnorm(i);
end
TS=TS';
figure
plotTS
RegionsLightArclight(:,3)=mean(TS(:,1000:(1000+50)),2)-mean(TS(:,(1000-50):1000),2);
RegionsOdorArclight(:,3)=mean(TS(:,7500:(7500+50)),2)-mean(TS(:,(7500-50):7500),2);
figure
plot(RegionsLightArclight,'+')
figure
plot(RegionsOdorArclight,'+')

figure
hold on
for i=1:3
plot(RegionsLightArclight(:,i)/mean(RegionsLightArclight(:,i)),'+')
end
figure
hold on
for i=1:3
plot(RegionsOdorArclight(:,i)/mean(RegionsOdorArclight(:,i)),'+')
end