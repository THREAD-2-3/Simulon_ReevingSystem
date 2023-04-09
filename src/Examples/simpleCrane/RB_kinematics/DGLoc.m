function DG=DGLoc(dq);

DG=zeros(3,4);

DG(1,1)=-dq(5);
DG(2,1)=-dq(6);
DG(3,1)=-dq(7);

DG(1,2)=dq(4);
DG(2,2)=-dq(7);
DG(3,2)=dq(6);

DG(1,3)=dq(7);
DG(2,3)=dq(4);
DG(3,3)=-dq(5);

DG(1,4)=-dq(6);
DG(2,4)=dq(5);
DG(3,4)=dq(4);

DG=2*DG;