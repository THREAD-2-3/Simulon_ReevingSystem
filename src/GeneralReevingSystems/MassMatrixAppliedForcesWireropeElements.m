function [M,Q] = MassMatrixAppliedForcesWireropeElements(p,dp,reevSys,ii)

ncwr = reevSys.ncwr; % Number of coordinates per element
nccr = reevSys.nccr; % Number of coordinates per element

M = zeros(ncwr*reevSys.nwr + nccr*reevSys.ncr,ncwr*reevSys.nwr + nccr*reevSys.ncr);
Q = zeros(ncwr*reevSys.nwr + nccr*reevSys.ncr,1);

for i = 1:(reevSys.nwr + reevSys.ncr)
    elemCoords = reevSys.elemCoords{i};
    elemCoords2 = elemCoords - reevSys.nrb*ones(1,length(elemCoords));
    
    if strcmp(reevSys.elementType(i),'wireRope')
        qe = p(elemCoords,1);
        dqe = dp(elemCoords,1);
        
        MyQ = MassMatrixAppliedForces(qe,dqe,reevSys.ro_wr(i));
        Me = MyQ(1:ncwr,1:ncwr);
        Qv = MyQ(1:ncwr,ncwr+1);
        Qgrav = MyQ(1:ncwr,ncwr+2);
        
        Qel = numericElasticForces(qe,dqe,reevSys.EA_wr(i),reevSys.EI_wr(i),reevSys.nu_wr(i),reevSys);
        
        M(elemCoords2,elemCoords2) = Me;
        Q(elemCoords2,1) = Qv + Qgrav + Qel;

%         if (ii==1500)&&(i==1)
%             [Qv Qgrav Qel]
%             [qe dqe]
%         end
    end
    
    if strcmp(reevSys.elementType(i),'cubicRope')
        qe = p(elemCoords,1);
        dqe = dp(elemCoords,1);
        
        MyQ = MassMatrixAppliedForces_cubicElement(qe,dqe,reevSys.ro_wr(i));
        Me = MyQ(1:nccr,1:nccr);
        Qv = MyQ(1:nccr,nccr+1);
        Qgrav = MyQ(1:nccr,nccr+2);
       
        %Qel = numericElasticForces_cubicElement_axialOnly(i,qe,reevSys.EA_wr(i),reevSys.EI_wr(i),reevSys);
        %Qel = numericElasticForces_cubicElement_e0(qe,reevSys.EA_wr(i),reevSys.EI_wr(i),reevSys);
        Qel = numericElasticForces_cubicElement_PVW(i,qe,reevSys.EA_wr(i),reevSys.EI_wr(i),reevSys);

        M(elemCoords2,elemCoords2) = Me;
        Q(elemCoords2,1) = Qv + Qgrav + Qel;
    end
end

% for i = 1:reevSys.nwr
%     qe = qwr(ncwr*(i-1)+1:ncwr*i,1);
%     dqe = dqwr(ncwr*(i-1)+1:ncwr*i,1);
%     
%     MyQ = MassMatrixAppliedForces(qe,dqe,reevSys.ro_wr(i));
%     Me = MyQ(1:ncwr,1:ncwr);
%     Qv = MyQ(1:ncwr,ncwr+1);
%     Qgrav = MyQ(1:ncwr,ncwr+2);
% 
%     Qel = numericElasticForces(qe,dqe,reevSys.EA_wr(i),reevSys.EI_wr(i),reevSys.nu_wr(i),reevSys);
% 
%     M(ncwr*(i-1)+1:ncwr*i,ncwr*(i-1)+1:ncwr*i) = Me;
%     Q(ncwr*(i-1)+1:ncwr*i,1) = Qv + Qgrav + Qel;
% end




