function G = GLoc(q)

G=zeros(3,4);

G(1,1)=-q(5);
G(2,1)=-q(6);
G(3,1)=-q(7);

G(1,2)=q(4);
G(2,2)=-q(7);
G(3,2)=q(6);

G(1,3)=q(7);
G(2,3)=q(4);
G(3,3)=-q(5);

G(1,4)=-q(6);
G(2,4)=q(5);
G(3,4)=q(4);

G=2*G;
