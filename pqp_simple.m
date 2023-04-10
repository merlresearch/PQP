function [x,values]=pqp(Q,h,x,iters,maxval)
% [x,values]=pqp(Q,h,x,iters,maxval) minimize f(x)==x'*Q*x/2-h'*x with x>=0
% INPUTS:
%  Q      -- symmetric positive semidefinite matrix 
%  h      -- vector
%  x      -- OPTIONAL initial guess
%  iters  -- OPTIONAL max # of iterations (default: 400)
%  maxval -- OPTIONAL maximum value(s) for x (box constraint)
% OUTPUTS:
%  x      -- final iterate of PQP multiplicative update
%  values -- OPTIONAL history of values of f(x)
%
% Copyright (C) 2008-2014,2023 Mitsubishi Electric Research Laboratories (MERL)
% SPDX-License-Identifier: AGPL-3.0-or-later

%%% Contents:
  %   About this algorithm/implementation
  %   Matlab/Octave code


%%% About this algorithm/implementation:

  % This is a baseline implementation of the Parallel Quadratic
  % Programming (PQP) algorithm, which in its simplest form is a
  % 1-line multiplicative fixpoint for the non-negative vector x>=0
  % that minimizes the quadratic form f(x)==x'*Q*x/2-h'*x, where Q is
  % symmetric semi-definite and h is any finite real vector.  The
  % fixpoint works by splitting f(x) into two a difference of two
  % quadratic forms, both of which are guaranteed to have strictly
  % non-negative gradients. The elementwise ratio of these gradients
  % is used to scale the solution estimate x, also elementwise; so all
  % elements of x can be updated in parallel, synchronously or
  % asynchronously, via very simple operations: two inner products,
  % one scalar division, and one scalar multiply.   This
  % multiplicative update is guaranteed to improve the value of f(x).  For
  % positive definite Q, repeated updates will converge at a linear
  % rate to the optimal x from any strictly positive initial guess.
  % For semidefinite Q, the rate of convergence will depend on the
  % splitting. Three of many possible splittings are implemented below
  % and can be chosen by uncommenting lines of code.  The simplest
  % split yields a 1-line algorithm:
  %
  % iterate: x = x.*(max(-Q,0)*x+max(h,0))./(max(Q,0)*x+max(-h,0));
  %
  % This code can also be used to solve general inequality-constrained
  % quadratic programs in their KKT dual form.  The baseline algorithm
  % can be quite fast and can solve large problems.  
  %
  % Applications, analyses, and accelerations can be found in:
  %
  % @inproceedings{Brand2011sep,
  % author = {Brand, M. and Chen, D.},
  % title = {Parallel Quadratic Programming for Image Processing},
  % booktitle = {IEEE International Conference on Image Processing (ICIP)},
  % year = 2011,
  % pages = {2261--2264},
  % month = sep,
  % doi = {10.1109/ICIP.2011.6116089},
  % url = {http://www.merl.com/publications/TR2011-064}
  % }
  %
  % @article{DiCairano2013jul,
  % author = {{Di Cairano}, S. and Brand, M. and Bortoff, S.A.},
  % title = {Projection-free Parallel Quadratic Programming 
  %          for Linear Model predictive Control},
  % journal = {International Journal of Control},
  % year = 2013,
  % month = jul,
  % url = {http://www.merl.com/publications/TR2013-059}
  % }
  %
  % @inproceedings{DiCairano2013dec2,
  % author = {{Di Cairano}, S. and Brand, M.},
  % title = {On a Multiplicative Update Dual Optimization Algorithm
  %          for Constrained Linear MPC},
  % booktitle = {IEEE Conference on Decision and Control (CDC)},
  % year = 2013,
  % month = dec,
  % url = {http://www.merl.com/publications/TR2013-108}
  % }
  %


%%% Matlab/Octave code:

%% INPUT CHECKING
if nargin<2, error('Need at least matrix Q and vector h'); end;
N=length(h); 
if [N,N]~=size(Q), error('Q and h are inconsistent sizes'); end
%% ALGORITHM PARAMETERS
epsilon=1e-6;			% a negligible value used to avoid 1/0 errors
threshold=N*epsilon;		% convergence threshold
%% SPLIT 
Qp=max(Q,0); hp=max(h,0);	% start with simplest possible split
%% ...CONTRACTION GUARANTEE: modify the split for semidefinite problems
% do nothing (fastest) or uncomment 1 of these 2 lines:
%Qp(1:N+1:end)=Qp(1:N+1:end)+sum(Qp-Q);  % (depends on Q) via diagonal dominance
%Qp=abs(Q);		   % (depends on Q and x)
Qp(1:N+1:end)=Qp(1:N+1:end)+epsilon; hp=hp*(1+epsilon); % improve conditioning
Qn=Qp-Q; hn=hp-h;					% complete split
%% INITIALIZATION
if (nargin<3)||(isempty(x)), x=(mean(abs(h))+abs(h))./diag(Q); end; % guess x
if (nargin<4)||(isempty(iters)), iters=1e4; end; % set default max iterations
if nargout>1, values(1)=0.5*(x'*Q*x)-h'*x; end; % user wants objective values
% ALGORITHM
for i=1:iters,
  x=max(x,epsilon); % OPTIONAL: bound x away from zero to avoid boundary traps
  x=x.*(Qn*x+hp)./(Qp*x+hn);	% PQP fixpoint (this is the algorithm)
  if nargin>=5, x=min(x,maxval); end; % OPTIONAL: force upper bound
  if nargout>1, values(i+1)=0.5*(x'*Q*x)-h'*x; end % OPTIONAL: record objective
  if mod(i,32)==0 && max(x.*abs(Q*x-h))<threshold, break; end % KKT converged?
end;
