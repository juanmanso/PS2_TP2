syms px py vx vy ax ay c11 c12 c21 c22 bx by T w

aa = (ax - bx) * T^2 / 2;
ab = (ay - by) * T^2 / 2;
ac = ((ax - bx) * T - (ay - by) * w * T^2 / 2);
ad = ((ay - by) * T - (ax - bx) * w * T^2 / 2);
s1 = 1 - w^2 * T^2 / 2;
s2 = w * T;

f1 = px + T * vx + aa * c11 + ab * c12;
f2 = py + T * vy + aa * c21 + ab * c22;
f3 = vx + ac * c11 + ad * c12;
f4 = vy + ac * c21 + ad * c22;
f5 = s1 * c11 + s2 * c12;
f6 = - s2 * c11 + s1 * c12;
f7 = s1 * c21 + s2 * c22;
f8 = - s2 * c21 + s1 * c22;
f9 = bx;
f10 = by;

f = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10];

Ad = jacobian(f, [px, py, vx, vy, c11, c12, c21, c22, bx, by])