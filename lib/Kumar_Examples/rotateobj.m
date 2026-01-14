function [xn, yn, zn] = rotateobj(x, y, z, rotmat)

for j = 1:size(x,1)
    for k = 1:size(x,2)
        vold = [x(j,k); y(j,k); z(j,k)];
        vnew = rotmat*vold;
        xn(j,k) = vnew(1); yn(j,k) = vnew(2); zn(j,k) = vnew(3);
    end
end