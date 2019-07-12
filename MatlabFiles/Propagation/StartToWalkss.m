syms PS(t) PB(t) VLP(t) SOG(t) AL(t) MB(t) LH(t)

a=1;
b=1;
c=1;
d=1;
e=1;
f=1;
g=0.2;
i=1;
j=1;
Walk=1;


ode1=diff(PS)==a*Walk;
ode2=diff(PB)==b*PS;
ode3=diff(VLP)==c*PS;
ode4=diff(SOG)==d*Walk;
ode5=diff(AL)==e*SOG;
ode6=diff(MB)==f*AL+g*SOG;
ode7=diff(LH)==i*AL+j*MB;

odes=[ode1; ode2; ode3; ode4; ode5; ode6; ode7];

cond1=PS(0)==0;
cond2=PB(0)==0;
cond3=VLP(0)==0;
cond4=SOG(0)==0;
cond5=AL(0)==0;
cond6=MB(0)==0;
cond7=LH(0)==0;

conds=[cond1;cond2;cond3;cond4;cond5;cond6;cond7];

S=dsolve(odes,conds);

figure
ezplot(S.PS)
hold on
ezplot(S.PB)
ezplot(S.VLP)
ezplot(S.SOG)
ezplot(S.AL)
ezplot(S.MB)
ezplot(S.LH)


