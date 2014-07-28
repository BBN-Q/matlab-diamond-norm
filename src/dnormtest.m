function dnormtest()
%DNORMTEST Test for diamond norm function
%
%   dnormtest() runs a series of minimal tests on the functionality of the diamond norm.
%
%
%  Copyright 2014 Raytheon BBN Technologies
%
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0                                  
%
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%

import qip.*
import qip.random.*
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

u2 = qip.random.unitary(2); uu2 = liou(u2,u2');
u4 = qip.random.unitary(4); uu4 = liou(u4,u4');

rr2 = @(c) choi_liou_involution(choi_liou_involution(c)*uu2);
rr4 = @(c) choi_liou_involution(choi_liou_involution(c)*uu4);

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

ct =[ dnorm(rr2(id)-rr2(id)),
      dnorm(rr2(id)-rr2(x)),
      dnorm(rr2(id)-rr2(y)),
      dnorm(rr2(id)-rr2(z)),
      dnorm(rr2(x)-rr2(y)),
      dnorm(rr2(x)-rr2(z)),
      dnorm(rr2(y)-rr2(x)),
      dnorm(rr2(y)-rr2(z)),
      dnorm(rr2(z)-rr2(x)),
      dnorm(rr2(z)-rr2(y)),
      dnorm(rr4(cnot)-rr4(id4)),
      dnorm(rr2(p1-p2)) ];

%t
%ct
%[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]'

%norm(t - [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]','inf')
%norm(ct - [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]','inf')

if (norm(t - [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]','inf') < 1e-7) & (norm(ct - [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sum(abs(r1-r2))/2]','inf') < 1e-7),
  disp('Test passed.');
else
  disp('Test failed.');
end



