function Fot = AnimateCrane(t,p,reevSys)

persistent hMS hDS hRB ha hb

[alpha,dalpha,ddalpha] = AzimutConstraint(t,reevSys);
[sa1,dsa1,ddsa1] = CoilingConstraint(t,reevSys);

A = eye(3);
A(1,1) = cos(alpha);
A(2,1) = sin(alpha);
A(1,2) = -sin(alpha);
A(2,2) = cos(alpha);

alf = 0:0.3:2*pi+0.1;

% Motor sheave

for i = 1:length(alf)
    rMS(:,i) = A*reevSys.R*[cos(alf(i)) 0.0 sin(alf(i))]';
end

% Deviating sheave

for i = 1:length(alf)
    rDS(:,i) = A*([reevSys.d 0 0]' + reevSys.R*[cos(alf(i)) 0 sin(alf(i))]');
end

% Rigid block

RRB = p(1:3,1);
ARB = rotMat(p(1:7,1));

uRB1 = [reevSys.b reevSys.b reevSys.b]';
uRB2 = [reevSys.b -reevSys.b reevSys.b]';
uRB3 = [-reevSys.b -reevSys.b reevSys.b]';
uRB4 = [-reevSys.b reevSys.b reevSys.b]';
uRB5 = [reevSys.b reevSys.b -reevSys.b]';
uRB6 = [reevSys.b -reevSys.b -reevSys.b]';
uRB7 = [-reevSys.b -reevSys.b -reevSys.b]';
uRB8 = [-reevSys.b reevSys.b -reevSys.b]';

rRB(:,1) = RRB + ARB*uRB1;
rRB(:,2) = RRB + ARB*uRB2;
rRB(:,3) = RRB + ARB*uRB3;
rRB(:,4) = RRB + ARB*uRB4;
rRB(:,5) = RRB + ARB*uRB1;
rRB(:,6) = RRB + ARB*uRB5;
rRB(:,7) = RRB + ARB*uRB6;
rRB(:,8) = RRB + ARB*uRB7;
rRB(:,9) = RRB + ARB*uRB8;
rRB(:,10) = RRB + ARB*uRB5;
rRB(:,11) = RRB + ARB*uRB8;


rRB(:,12) = RRB + ARB*uRB4;
rRB(:,13) = RRB + ARB*uRB3;
rRB(:,14) = RRB + ARB*uRB7;
rRB(:,15) = RRB + ARB*uRB3;
rRB(:,16) = RRB + ARB*uRB2;
rRB(:,17) = RRB + ARB*uRB6;

% Wire rope 'a'

qe = p(8:24,1);
s1 = qe(16); s2 = qe(17);
qa = qe(1:6,1);
s = s1:(s2-s1)/10:s2;

r1 = qa(1:3,1);
r2 = qa(4:6,1);
qm = qe(7:12,1);
ie = (r2-r1)/sqrt((r2-r1)'*(r2-r1));
ke = cross(ie,[0 1 0]');
je = cross(ke,ie);
Ae = [ie je ke];

for i = 1:length(s)
    shi = (2*s(i)-s1-s2)/(s2-s1);
    N1 = (1-shi)/2;
    N2 = (1+shi)/2;
    N = [N1*eye(3) N2*eye(3)];
    ra_3d = N*qa;
    
    S1 = sin(pi*(s(i)-s1)/(s2-s1));
    S2 = sin(2*pi*(s(i)-s1)/(s2-s1));
    S3 = sin(3*pi*(s(i)-s1)/(s2-s1));
    
    S(2,1) = S1; S(2,2) = S2; S(2,3) = S3;
    S(3,3+1) = S1; S(3,3+2) = S2; S(3,3+3) = S3;
    ut_3d = S*qm;
    
    ra(:,i) = ra_3d + ut_3d;
end

% Wire rope 'b'

qe = p(25:41,1);
s1 = qe(16); s2 = qe(17);
qa = qe(1:6,1);
s = s1:(s2-s1)/10:s2;

r1 = qa(1:3,1);
r2 = qa(4:6,1);
qm = qe(7:12,1);
ie = (r2-r1)/sqrt((r2-r1)'*(r2-r1));
ke = cross(ie,[0 1 0]');
je = cross(ke,ie);
Ae = [ie je ke];

for i = 1:length(s)
    shi = (2*s(i)-s1-s2)/(s2-s1);
    N1 = (1-shi)/2;
    N2 = (1+shi)/2;
    N = [N1*eye(3) N2*eye(3)];
    rb_3d = N*qa;
    
    S1 = sin(pi*(s(i)-s1)/(s2-s1));
    S2 = sin(2*pi*(s(i)-s1)/(s2-s1));
    S3 = sin(3*pi*(s(i)-s1)/(s2-s1));
    
    S(2,1) = S1; S(2,2) = S2; S(2,3) = S3;
    S(3,3+1) = S1; S(3,3+2) = S2; S(3,3+3) = S3;
    ut_3d = S*qm;
    
    rb(:,i) = rb_3d + ut_3d;
end

if t==0
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[0.1*scrsz(4) 0.1*scrsz(4) 0.3*scrsz(3) 0.8*scrsz(4)])
    
    hMS = plot3(rMS(1,:),rMS(2,:),rMS(3,:),'b','XDataSource','rMS(1,:)','YDataSource','rMS(2,:)','ZDataSource','rMS(3,:)','LineWidth',2);
    hold on;
    hDS = plot3(rDS(1,:),rDS(2,:),rDS(3,:),'b','XDataSource','rDS(1,:)','YDataSource','rDS(2,:)','ZDataSource','rDS(3,:)','LineWidth',2);
    hRB = plot3(rRB(1,:),rRB(2,:),rRB(3,:),'g','XDataSource','rRB(1,:)','YDataSource','rRB(2,:)','ZDataSource','rRB(3,:)','LineWidth',2);
    ha = plot3(ra(1,:),ra(2,:),ra(3,:),'k','XDataSource','ra(1,:)','YDataSource','ra(2,:)','ZDataSource','ra(3,:)','LineWidth',2);
    hb = plot3(rb(1,:),rb(2,:),rb(3,:),'k','XDataSource','rb(1,:)','YDataSource','rb(2,:)','ZDataSource','rb(3,:)','LineWidth',2);
    
    axis equal;
    axis([-1, 5,-1, 5, -10, 1]);
    %axis([-0.5, 3.5, -0.1, 0.2]);
    %axis([2, 4, -9, -6]);
    %xlabel('Distance along the track (m)')
else
    axis equal;
    axis([-1, 5,-1, 5, -10, 1]);
    %axis([-0.5, 3.5, -0.1, 0.2]);
    %axis([2, 4, -9, -6]);
    
    refreshdata(hMS,'caller')
    refreshdata(hDS,'caller')
    refreshdata(hRB,'caller')
    refreshdata(ha,'caller')
    refreshdata(hb,'caller')
    drawnow;
end

Fot = getframe;