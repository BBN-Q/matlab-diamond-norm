import qip.*
import qip.open_systems.*

id = choi_liou_involution(liou(pauli(0),pauli(0)));
x  = choi_liou_involution(liou(pauli(1),pauli(1)));
y  = choi_liou_involution(liou(pauli(2),pauli(2)));
z  = choi_liou_involution(liou(pauli(3),pauli(3)));

r1 = rand(4); r1 = r1/sum(r1);
r2 = rand(4); r2 = r2/sum(r2);

p1 = r1(1)*id+r1(2)*x+r1(3)*y+r1(4)*z;
p2 = r2(1)*id+r2(2)*x+r2(3)*y+r2(4)*z;

cnot_u = [eye(2) zeros(2); zeros(2) pauli(1) ];
cnot = choi_liou_involution(liou(cnot_u, cnot_u));
id4  = choi_liou_involution(liou(eye(4), eye(4)));

t = [ dnorm(id-id),
      dnorm(id-x),
      dnorm(id-y),
      dnorm(id-z),
      dnorm(x-y),
      dnorm(x-z),
      dnorm(y-x),
      dnorm(y-z),
      dnorm(z-x),
      dnorm(z-y),
      dnorm(cnot-id4),
      dnorm(p1-p2) ];

if norm(t - [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]','inf') < 1e-11,
  disp('Test passed.');
else
  disp('Test failed.');
end


