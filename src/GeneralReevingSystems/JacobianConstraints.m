function Cp = JacobianConstraints(p,reevSys)

Cp = JacobNodalCoordinatesConstraints(p,reevSys);

Cp_RB = RigidBodyJacobianConstraints(p,reevSys);

Cp = [Cp_RB' Cp']';
