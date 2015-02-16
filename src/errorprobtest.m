function pass=errorprobtest(N)
import qip.*
import qip.random.*
import qip.open_systems.*

cpauli = @(n) qip.open_systems.choi_liou_involution(qip.open_systems.liou(qip.pauli(n),qip.pauli(n)));
pass = true;

disp('Testing Pauli channels ...')
% try a bunch of random pauli channels
for jj=1:N
  p = rand(4,1); p = p/sum(p);
  c = zeros(4,4);
  for ii=1:4,
    c = c + p(ii)*cpauli(ii-1);
  end
  if abs(errorprob(c)-sum(p(2:4)))>1e-6
    disp('Test failed')
    pass = false;
    break
  end
end

disp('Testing random CP maps ...')
% try a bunch of random CP maps
for jj=1:N
  p = rand(2,1); p = p/sum(p);
  c = p(1)*cpauli(0)+p(2)*choi_liou_involution(choi2liou(choi_matrix(2)));
  if errorprob(c)-p(2)>0
    disp('Test failed')
    disp(p)
    disp(c)
    disp(errorprob(c))
    pass = false;
    break
  end
end

disp('Testing relaxation ...')
% relaxation is non stochastic
for jj=1:N
  p = 1000*rand();
  q = rand(2,1);q=q/sum(q);
  c = q(1)*cpauli(0)+q(2)*choi_liou_involution(expm(1/p*dissipator([0 1; 0 0])));
  if abs(q(2)-errorprob(c))>1e-3 % require 'precision' in errorprob to be set to 'high'
    disp('Test failed')
    disp(p)
    disp(q)
    disp(c)
    disp(errorprob(c))
    pass = false;
    break
  end
end

if pass
  disp('All tests passed')
end
return
