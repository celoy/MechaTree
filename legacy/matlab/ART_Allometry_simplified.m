Rn=4;
D=2;
k=(1:1:10);

S=Rn.^(k*0.1);
l= (Rn.^(k/D)-1) / (Rn.^(1/D)-1);
diam=S.^(1/2).*(1+0.2*k).^(1/3).*Rn.^(k/2);

diam=( l.*( 1+l/l(size(l,2)) ) ).^(1/3);

loglog(diam,l,'o-')

%%
l=(1:100);
d=l.^2.*sqrt(log(l));

loglog(d,l)