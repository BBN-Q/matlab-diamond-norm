function result = errordnorm( J )
%ERRORDNORM Diamond norm of error
%
%   errordnorm(J) returns the diamond form of the difference between a
%            Choi matrix J and the Choi matrix corresponding to the identity channel
%            of the appropriate dimension.
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
  
dim = round(sqrt(size(J,1)));

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

delta = J - choi_id;

result = dnorm(delta);
end

