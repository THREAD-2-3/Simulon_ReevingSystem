function C = Constraints(p,reevSys)

C_WR = NodalCoordinatesConstraints(p,reevSys);
C_RB = RigidBodyConstraints(p,reevSys);

C = [C_RB' C_WR']';