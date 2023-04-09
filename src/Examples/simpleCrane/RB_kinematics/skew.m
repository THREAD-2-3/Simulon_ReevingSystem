function vsk = skew(v)

vsk=zeros(3,3);

vsk(1,2)=-v(3);
vsk(1,3)=v(2);
vsk(2,3)=-v(1);

vsk(2,1)=-vsk(1,2);
vsk(3,1)=-vsk(1,3);
vsk(3,2)=-vsk(2,3);