function result = errorprob( J )
%PERROR Error probability
%
%   errorprob(J) returns the probability of a stochastic error in the superoperator
%     represented by the Choi matrix J. In essence, it looks for the smallest p such 
%     that 
%
%          J - (1-p)*I >= 0, 
%
%     where I is the Choi matrix for the identity channel.
%
%   The normalization of J must be such that partial prace of the first
%   half of J always yields the identity for trace preserving
%   superoperators.
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
  
dimsq = size(J,1);
dim = round(sqrt(dimsq));
Jh = (J'+J)/2;

persistent choi_id;
persistent old_dim;

if size(choi_id) == [0 0] | dim ~= old_dim,
  old_dim = dim;
  max_ent = zeros(dim^2,1);
  unit = zeros(dim,1); unit(1) = 1;
  for ii=1:dim
    max_ent = max_ent + kron(unit,unit);
    unit = circshift(unit,1);
  end
  choi_id = max_ent*max_ent';
end

cvx_begin sdp quiet
  cvx_precision default % replace 'default' with 'high' for 1e-10 precision, but slower calculation

  variable p;

  maximize( p );
  subject to
    p >= 0;
    p <= 1;
    Jh - p*choi_id >= 0;

cvx_end

result = 1-cvx_optval;

end

