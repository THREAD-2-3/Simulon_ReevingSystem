function w = Omega(q,dq)

G = Gglob(q);

w = G*dq(4:7,1);