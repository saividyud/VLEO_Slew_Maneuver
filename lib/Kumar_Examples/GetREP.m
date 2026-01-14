function R = GetREP(be, Ph)
cp = cos(Ph); 
sp = sin(Ph);
lam = 1 - cp;

b12 = be(1)*be(2)*lam;
b13 = be(1)*be(3)*lam;
b23 = be(2)*be(3)*lam;

R = [(be(1)^2*lam + cp), (b12 + be(3)*sp), (b13 - be(2)*sp);
    (b12 - be(3)*sp), (be(2)^2*lam  + cp), (b23 + be(1)*sp);
    (b13 + be(2)*sp), (b23 - be(1)*sp), (be(3)^2*lam + cp)];
