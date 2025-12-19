function result = band_Lanczos(A,R,L,nmax,sflag,tol,n,result)
%
%  Band Lanczos method
% 
%  This implementation is essentially Algorithm 5.1 in
% 
%  Roland W. Freund, Model reduction methods based on Krylov subspaces, 
%  Acta Numerica, 12 (2003), pp. 267-319.
%
% -----------------------------------------------------------------------------
%
%  Notes:  1) The matrices A, R, and L are allowed to be complex.
%
%          2) The biorthogonalizations performed in the algorithm are with
%             respect to the bilinear form
%  
%                    (w,v) = w^T v (= w.' * v),
%
%             where w^T (= w.') denotes the transpose of w.
%
%          3) This algorithm explicitly enforces only the minimum possible 
%             amount of biorthogonalization.  Moreover, only the entries below
%             (and for T on) the diagonal of the Lanczos matrices T and Tt are 
%             computed directly via inner products and norms of candidate
%             vectors and previously deflated vectors;  the other entries are
%             obtained via the relation
% 
%                    D * T = Tt.' * D
%
%             between T and the transpose Tt.' of Tt.
% 
%          4) This algorithm has no built-in look-ahead procedure to continue
%             in the case of a breakdown or near-breakdown, and instead, the
%             algorithm simply stops in this case. 
%
% -----------------------------------------------------------------------------
%
%  Usage:  result = band_Lanczos(A,R,L,nmax)
%          result = band_Lanczos(A,R,L,nmax,sflag)
%          result = band_Lanczos(A,R,L,nmax,sflag,tol)
%          result = band_Lanczos(A,R,L,nmax,sflag,tol,n,result)
% 
%          where A is a matrix, or,
%
%          result = band_Lanczos(@(x,tfl) afun(x,tfl,...),R,L,nmax)
%          result = band_Lanczos(@(x,tfl) afun(x,tfl,...),R,L,nmax,sflag)
%          result = band_Lanczos(@(x,tfl) afun(x,tfl,...),R,L,nmax,sflag,tol)
%          result = band_Lanczos(@(x,tfl) afun(x,tfl,...),R,L,nmax,sflag, ...
%                                                                tol,n,result)
%
%          where "afun" is a function such that 
% 
%          y = afun(x,0,...)
%
%          computes the matrix-vector product y = A x (= A * x) with A and
%
%          y = afun(x,1,...)
%
%          computes the matrix-vector product y = A^T x (= A.' * x) with the
%          transpose A^T  of A 
%
% -----------------------------------------------------------------------------
%         
%  Required inputs:  A = a square matrix or an (anonymous) function that 
%                        computes matrix-vector products with A and A^T
%                    R = matrix the m columns of which are the right starting
%                        vectors
%                    L = matrix the p columns of which are the left starting
%                        vectors
%                 nmax = maximum number of pairs of right and left Lanczos
%                        vectors to be generated
%
%                 It is assumed that A is a square matrix and that A, R, and L
%                 have the same number of rows;  these assumptions are checked
%                 when the input A is a matrix, but not when A is a function. 
%
% -----------------------------------------------------------------------------
%
%  Optional inputs:  sflag = an integer that controls the storage of the right
%                            and left Lanczos vectors:  sflag = 0 means that 
%                            only the vectors needed for the
%			     biorthogonalization are stored;  a nonzero value 
%                            of sflag means that all vectors are stored
%                      tol = a structure that contains tolerances and
%			     parameters for the deflation procedure and the
%             		     breakdown check
%                        n = a nonnegative integer;  n > 0 means that a
%                            previous call to the function "band_Lanczos" has
%                            generated n pairs of right and left Lanczos
%                            vectors and that the current call to 
%                            "band_Lanczos" resumes the iteration at step n+1;
%                            n = 0 means that the band Lanczos method is
%                            started from scratch
%                   result = the output structure of the previous call to
%                            "band_Lanczos" if n > 0;  if n = 0, this input is
%                            ignored
% 
%                 If sflag is not provided as an input, sflag = 0 is used.
%
%                 If "tol" is provided as an input, it needs to contain 
%
%                    tol.defl_tol = unscaled deflation tolerance
%
%                    tol.brk_tol = tolerance for breakdown check
%
%                 and values of the following two flags:
%
%                    tol.defl_flag = 0  (use unscaled deflation tolerance)
%                                    1  (use scaled deflation tolerance)
%                                    2  (use scaled deflation tolerance only
%                                            for the starting blocks R and L)
%                                    3  (use scaled deflation tolerance except
%                                              for the starting blocks R and L)
%
%                    tol.normA_flag = 1  (an estimate for the norm of A is
%                                           generated within the algorithm)
%                                     2  (an estimate for the norm of A is
%                                                    provided as tol.normA)
% 
%                 If tol.normA_flag = 2, then an estimate for the norm of A 
%                 needs to be provided as 
%
%                    tol.normA
%
%                 If "tol" is not provided as an input, the following default
%                 values are used:
%
%                    tol.defl_tol = sqrt(eps)  (eps = machine precision)
%                    tol.brk_tol = 100 * eps;
%                    tol.defl_flag = 1
%                    tol.normA_flag = 1
%
%                 If n is not provided as an input, n = 0 is used and the
%                 optional input "result" is ignored.
%							 a
%                 If n > 0, the input "result" needs to be provided.  It is
%                 assumed, but not checked, that "result" is the output 
%                 structure of a previous call to the function "band_Lanczos"
%                 (applied to the same matrices A, R, and L).  In this case, 
%                 the inputs "sflag" and "tol" are ignored and
%                  
%                    tol = result.tol  
%
%                 and
%
%                    sflag = result.sflag
%
%                 are used instead.
%
% -----------------------------------------------------------------------------
%
%  On return, "result" is a structure the fields of which include the following
%  quantities:
%
%         result.n = number of pairs of right and left Lanczos vectors that
%                    were generated
%
%         result.V = matrix V the columns of which are the stored right Lanczos
%                    vectors;  if sflag = 0, V has at most p+1 columns and
%		     these columns are the right Lanczos vectors that are 
%                    needed for biorthogonalization;  if sflag is nonzero, V
%                    contains the first n right Lanczos vectors
%
%         result.W = matrix W the columns of which are the stored left Lanczos 
%                    vectors;  if sflag = 0, W has at most m+1 columns and
%		     these columns are the left Lanczos vectors that are needed
%		     for biorthogonalization;  if sflag is nonzero, W contains
%                    the first n left Lanczos vectors
%       
%   result.Vh_defl = matrix Vh_defl the columns of which are the candidate
%                    vectors for the next mc right Lanczos vectors and the m-mc
%		     deflated right vectors 
%      
%   result.Wh_defl = matrix Wh_defl the columns of which are the candidate
%                    vectors for the next pc left Lanczos vectors and the p-pc
%                    deflated left vectors 
%      
%         result.D = vector D the entries of which are the products
%
%                             (W(:,k)).' * V(:,k),  k = 1, 2, ..., n
%
%         result.T = the n x n matrix T that represents the oblique projection
%                    of the matrix A onto the subspace spanned by the columns 
%                    of V, and orthogonally to the subspace spanned by the
%                    columns of W;  A, V, W, and T are connected via the
%                    relation
%
%                             W.' * A * V = (W.' * V) * T
%
%        result.Tt = the n x n matrix Tt that represents the oblique projection
%                    of the matrix A.' onto the subspace spanned by the columns
%                    the columns of W, and orthogonally to the subspace spanned
%                    by the columns of V;  A, V, W, and Tt are connected via  
%                    the relation 
%
%                             V.' * A.' * W = (V.' * W) * Tt
%
%       result.rho = the matrix rho that contains the coefficients used to turn
%                    the right starting vectors (in R) into the first right 
%                    Lanczos vectors;  R, V, W, and rho are connected via the
%                    relation 
% 
%                             W.' * R = W.' * V * rho
%
%       result.eta = the matrix eta that contains the coefficients used to turn
%                    the left starting vectors (in L) into the first left
%                    Lanczos vectors;  L, V, W, and eta are connected via the
%                    relation 
%
%                             V.' * L = V.' * W * eta
% 
% -----------------------------------------------------------------------------
%
%  This routine can be run in incremental fashion. 
%
%  Example 1:   
%
%      result = band_Lanczos(A,R,L,nmax0,sflag,tol)
%    
%      n = result.n					 
%
%      result = band_Lanczos(A,R,L,nmax,sflag,tol,n,result)
%
%      The first call to the function "band_Lanczos" runs the band Lanczos 
%      method from scratch and generates n pairs of right and left Lanczos 
%      vectors.
%  
%      The second call to the function of "band_Lanczos" resumes the iteration 
%      at step n+1.
% 
%  Example 2 (Band Lanczos method, run one step at a time):
%
%      result = band_Lanczos(A,R,L,1,sflag,tol,0,[])
% 
%      for n = 1 : nmax - 1,
% 
%         result = band_Lanczos(A,R,L,n+1,[],[],result.n,result)
%  
%      end
% 
%      This will run the band Lanczos method for nmax steps.
%
% -----------------------------------------------------------------------------
%
%  BANDITS: a Matlab Package of Band Krylov Subspace Iterations
%
%  Copyright (c) 2018-2019 Roland W. Freund
%  See LICENSE.txt for license
%
% -----------------------------------------------------------------------------

%  =======================
%  Begin of initialization
%  =======================
%
if nargin < 4,
   error('** Not enough input arguments! **')
end
%
if (nargin < 5) || isempty(sflag),
   sflag = 0;
else
   if (rem(sflag,1) ~= 0),
      error('** sflag needs to be an integer **')
   end  
end  
%
if (nargin < 6) || isempty(tol),
   tol.defl_flag = 1;
   tol.defl_tol = sqrt(eps);
   tol.normA_flag = 1;
   tol.brk_tol = 100 * eps;
end
%
if (nargin < 7) || isempty(n),
   n = 0;
else
   if (n < 0) || (rem(n,1) ~= 0),
      error('** n needs to be a nonnegative integer **')
   end  
end
%
if nmax <= n,
   error('** nmax is not large enough;  we need to have nmax > n **')
end  
%
[N,m] = size(R);
[Nt,p] = size(L);
%
if N ~= Nt,
   error('** R and L need to have the same number of rows **')
end   
%
if isfloat(A) == 1,
   mvec_A = 1;
   [Nt1,Nt2] = size(A);
   if Nt1 ~= Nt2,
      error('** The matrix A is not square **')
   end
   if Nt1 ~= N,
      error('** The matrices A, R, L need to have the same number of rows **')
   end
else
   mvec_A = 0;
end
%
if n > 0,
   if (nargin < 8) || isempty(result),
      error('** n > 0, but there is no input "result" **')
   end   
   V = result.V;
   W = result.W;
   Vh_defl = result.Vh_defl;
   Wh_defl = result.Wh_defl;  
   D = result.D;
   T = result.T;
   Tt = result.Tt;  
   rho = result.rho;
   eta = result.eta;
   mc = result.mc;
   pc = result.pc;
   Iv = result.Iv;
   Iw = result.Iw;
   n_check = result.n;
   sflag = result.sflag;
   tol = result.tol;
   brk_flag = result.brk_flag;
   exh_flag = result.exh_flag;
%
   if brk_flag == 1,
      fprintf(' \n')
      disp('**-------------------------------------------**')
      disp('** Previous run ended with breakdown;        **')
      disp('** look-ahead is needed in order to continue **')    
      disp('**-------------------------------------------**')
      fprintf(' \n')    
      return
   end
%  
   if exh_flag > 0,
      fprintf(' \n')
      disp('**-----------------------------------------------**')    
      disp('** Previous run ended due to an                  **')    
      disp('** exhausted right or left block Krylov subspace **')    
      disp('**-----------------------------------------------**')
      fprintf(' \n') 
      return
   end
%
   if n ~= n_check,
      error('** n does not match the value of n in result **')
   end
%  
   n1 = n + 1;
%
else
   V = zeros(N,0);
   W = zeros(N,0);  
   Vh_defl(:,1:m) = R;
   Wh_defl(:,1:p) = L;
   D = [];
   T = [];
   Tt = [];
   rho = [];
   eta = [];
   mc = m;
   Iv.v = [];
   if sflag == 0,
      Iv.av = [1:p+1];
      Iw.av = [1:m+1];
   end
   Iv.ph = [1:m];
   Iv.I = [];
   Iv.pd = [];  
   Iv.nd = 0;
   pc = p;
   Iw.v = [];
   Iw.ph = [1:p];
   Iw.I = [];
   Iw.pd = [];  
   Iw.nd = 0;
   n1 = 1;
   brk_flag = 0;
   exh_flag = 0;
end
%
result.m = m;
result.p = p;
%
%  Extract and check tolerances, flags, and norm estimate for deflation 
%  and breakdown checks
%
[defl_tol,defl_flag,normA_flag,normA,brk_tol] = check_tolerances(tol,n);
%
%  ============================================
%  End of initialization and begin of iteration
%  ============================================
%
for n = n1 : nmax,
%
%  ==================================================
%  Construct n-th pair of Lanczos vectors v_n and w_n
%  ==================================================
%
   foundvn = 0;
   foundwn = 0; 
%   
%  If necessary, deflate v or w vector
%
   while (foundvn == 0) | (foundwn == 0),
%
      if mc >= pc,
%         
%        If mc >= pc, we first check for deflation of v vectors
%
         if foundvn == 0,
%         
            [mc,foundvn,Vh_defl,Iv,normv] = deflation(1,n,m,mc,foundvn, ...
                               Vh_defl,R,Iv,defl_flag,defl_tol,normA);
%      
         end
%       
         if foundwn == 0,
%
            [pc,foundwn,Wh_defl,Iw,normw] = deflation(-1,n,p,pc,foundwn, ...
                                Wh_defl,L,Iw,defl_flag,defl_tol,normA);
%               
         end
%        
      else
%
%        In this case pc > mc, and we first check for deflation of w vectors
%
         if foundwn == 0,
%
            [pc,foundwn,Wh_defl,Iw,normw] = deflation(-1,n,p,pc,foundwn, ...
                                Wh_defl,L,Iw,defl_flag,defl_tol,normA);
%
         end
%         
         if foundvn == 0,
%         
            [mc,foundvn,Vh_defl,Iv,normv] = deflation(1,n,m,mc,foundvn, ...
                               Vh_defl,R,Iv,defl_flag,defl_tol,normA);
%       
         end
%     
%     End of:  if mc >= pc,
%
      end
%       
%     Check if block Krylov subspace is exhausted
%
      if (mc == 0) | (pc == 0),
%             
         disp(['  Number of Lanczos steps performed: ' num2str(n-1)])
%
         tol.normA = normA;
         result = save_result(1,V,Vh_defl,mc,Iv,T,rho,tol,D, ...
                                            W,Wh_defl,pc,Iw,Tt,eta);
%
         result.n = n - 1;
         result.sflag = sflag;      
         result.brk_flag = brk_flag;
%         
         if mc == 0,
            disp('**---------------------------------------------**')
            disp('** There are no more right Krylov vectors,     **')
            disp('** and so the algorithm has to terminate: STOP **')
            disp('**---------------------------------------------**')
            result.exh_flag = 1;
         end
%             
         if pc == 0,
            disp('**---------------------------------------------**')
            disp('** There are no more left Krylov vectors,      **')
            disp('** and so the algorithm has to terminate: STOP **')
            disp('**---------------------------------------------**')
            result.exh_flag = 2;             
         end
%             
         return
%      
      end
%
%  End of:  while (foundvn == 0) | (foundwn == 0),  
%
   end
%
%  Normalize v_n and w_n
%
   tmpv = Vh_defl(:,Iv.ph(1)) / normv;
   tmpw = Wh_defl(:,Iw.ph(1)) / normw;
%  
%  Compute delta_n and check for breakdown
%
   delta = tmpw.' * tmpv;
%
   if abs(delta) <= brk_tol,
      disp('**------------------------------- **')
      disp('** Breakdown --- look-ahead would **')
      disp('** be needed to continue: STOP    **')
      disp('**--------------------------------**')
      disp(['  Number of Lanczos steps performed: ' num2str(n-1)])
      disp(['  Tolerance for breakdown check:     ' num2str(brk_tol)])
      disp(['  Absolute value of delta = w^T * v: ' num2str(abs(delta))])
      disp('**---------------------------------------------**'), fprintf(' \n')
%
      tol.normA = normA;
      result = save_result(1,V,Vh_defl,mc,Iv,T,rho,tol,D, ...
                                         W,Wh_defl,pc,Iw,Tt,eta);
%
      result.n = n - 1;
      result.sflag = sflag;      
      result.brk_flag = 1;
      result.exh_flag = 0;      
      return      
   end
%  
   D(n,1) = delta;
%    
%  Put v_n and w_n into their slots in V and W, respectively
%
   if sflag == 0,
      nvi = min(Iv.av);
      Iv.av = setdiff(Iv.av,nvi);
      nwi = min(Iw.av);
      Iw.av = setdiff(Iw.av,nwi);
   else
      nvi = n;
      nwi = n;
   end
%  
%  Make sure rho and eta have n rows
%
   rho(n,1) = 0;
   eta(n,1) = 0;    
%
   V(:,nvi) = tmpv;
   Iv.v(n) = nvi;
   if n > mc,
      T(n,n-mc) = normv;
   else
      rho(n,n-mc+m) = normv;
   end
%
   W(:,nwi) = tmpw;
   Iw.v(n) = nwi;
   if n > pc,
      Tt(n,n-pc) = normw;
   else
      eta(n,n-pc+p) = normw;
   end
%
%  Biorthogonalize the right candidate vectors against w_n
%
   ivph1 = Iv.ph(1);
   Itmp = Iv.ph(2:mc);
   Iv.ph(1:mc-1) = Itmp;
   Iv.ph(mc) = ivph1;
%
   tmp = ((W(:,nwi)).' * Vh_defl(:,Itmp)) / delta;
   Vh_defl(:,Itmp) = Vh_defl(:,Itmp) - V(:,nvi) * tmp; 
%
   Ktmp = find([1:mc-1] > mc-n);
   T(n,Ktmp-mc+n) = tmp(Ktmp);
%  
   Ktmp = find([1:mc-1] <= mc-n);
   rho(n,Ktmp-mc+n+m) = tmp(Ktmp);  
%
%  Biorthogonalize the left candidate vectors against v_n
%
   iwph1 = Iw.ph(1);
   Itmp = Iw.ph(2:pc);
   Iw.ph(1:pc-1) = Itmp;
   Iw.ph(pc) = iwph1;
%
   tmp = ((V(:,nvi)).' * Wh_defl(:,Itmp)) / delta;
   Wh_defl(:,Itmp) = Wh_defl(:,Itmp) - W(:,nwi) * tmp; 
%
   Ktmp = find([1:pc-1] > pc-n);
   Tt(n,Ktmp-pc+n) = tmp(Ktmp);
%  
   Ktmp = find([1:pc-1] <= pc-n);
   eta(n,Ktmp-pc+n+p) = tmp(Ktmp);  
%
%  Advance right block Krylov subspace by computing 
%  tmpv = A * V(:,nvi)) (= A v_n)
%
   if mvec_A == 1,
      tmpv = A * V(:,nvi);
   else
      tmpv = feval(A,V(:,nvi),0);
   end
%
   if normA_flag == 1,
      normA = max([normA, norm(tmpv,2)]);
   end
%
%  Biorthogonalize tmpv against w vectors
%    
   nd = Iw.nd;
   tmp = (V(:,nvi)).' * Wh_defl(:,Iw.pd(1:nd));
%    
   Itmp = Iw.I(1:nd);
   Ktmp = find(Itmp > 0);
   IKtmp = Itmp(Ktmp);
   Tt(n,IKtmp) = tmp(Ktmp) / delta;
   T(IKtmp,n) = (tmp(Ktmp)).' ./ D(IKtmp,1);
   tmpv = tmpv - V(:,Iv.v(IKtmp)) * T(IKtmp,n);
%
   Ktmp = find(Itmp <= 0);
   eta(n,Itmp(Ktmp)+p) = tmp(Ktmp) / delta;
%
   Ktmp = max(1,n-pc):n-1;
   T(Ktmp,n) = ((Tt(n,Ktmp)).' * delta) ./ D(Ktmp,1);
   tmpv = tmpv - V(:,Iv.v(Ktmp)) * T(Ktmp,n);
%
   T(n,n) = ((W(:,nwi)).' * tmpv) / delta;
   Vh_defl(:,Iv.ph(mc)) = tmpv - V(:,nvi) * T(n,n);
%
%  Advance left block Krylov subspace) by computing 
%  tmpw = A.' * W(:,nwi)) (=  A^T w_n)
%
   if mvec_A == 1,
      tmpw = A.' * W(:,nwi);
   else
      tmpw = feval(A,W(:,nwi),1);
   end
%
   if normA_flag == 1,
      normA = max([normA, norm(tmpw,2)]);
   end
%
%  Biorthogonalize tmpw against v vectors
%    
   nd = Iv.nd;
%    
   tmp = (W(:,nwi)).' * Vh_defl(:,Iv.pd(1:nd));
%    
   Itmp = Iv.I(1:nd);
   Ktmp = find(Itmp > 0);
   IKtmp = Itmp(Ktmp);
   T(n,IKtmp) = tmp(Ktmp) / delta;
   Tt(IKtmp,n) = (tmp(Ktmp)).' ./ D(IKtmp,1);
   tmpw = tmpw - W(:,Iw.v(IKtmp)) * Tt(IKtmp,n);
%
   Ktmp = find(Itmp <= 0);
   rho(n,Itmp(Ktmp)+m) = tmp(Ktmp) / delta;
%
   Ktmp = max(1,n-mc):n-1;
   Tt(Ktmp,n) = ((T(n,Ktmp)).' * delta) ./ D(Ktmp,1);
   tmpw = tmpw - W(:,Iw.v(Ktmp)) * Tt(Ktmp,n);
%   
   Tt(n,n) = T(n,n);
   Wh_defl(:,Iw.ph(pc)) = tmpw - W(:,nwi) * Tt(n,n);
%
%  Unless all Lanczos vectors are stored, adjust available slots
%  for right and left Lanczos vectors in V and W, respectively
%
   if sflag == 0,
      if n > pc,
         if (ismember(n-pc,Iw.I) == 0),
            Iv.av = union(Iv.av,Iv.v(n-pc));
         end
      end
%
      if n > mc,
         if (ismember(n-mc,Iv.I) == 0),
            Iw.av = union(Iw.av,Iw.v(n-mc));
         end    
      end
   end
%  
%  End of:  for n = n1 : nmax,
% 
end
%
%  ================
%  End of iteration
%  ================
%
tol.normA = normA;
result = save_result(1,V,Vh_defl,mc,Iv,T,rho,tol,D,W,Wh_defl,pc,Iw,Tt,eta);
%
result.n = n;
result.sflag = sflag;      
result.brk_flag = 0;
result.exh_flag = 0;

