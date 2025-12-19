function [mc,fv,Vh,I,nv] = deflation(vw_f,n,m,mc,fv,Vh,R,I,d_f,d_t,nA)
%
%  This function is an auxiliary routine that is used to check for deflation
%  and if necessary, to perform the deflation.  It is called in the functions
%  "band_Arnoldi", "band_Lanczos", "sym_band_Lanczos", "Herm_band_Lanczos",
%  "Jsym_band_Lanczos", and "JHerm_band_Lanczos".
%
%  The input flag "vw_f" is used to produce correct display messages:
%
%  vw_f = 0   means that there is only one type of Krylov vectors (as in 
%             "band_Arnoldi", "Herm_band_Lanczos", "sym_band_Lanczos",
%             "Jsym_band_Lanczos", "JHerm_band_Lanczos")
%
%  vw_f = 1   means that right Krylov vectors are deflated
%
%  vw_f = -1  means that left Krylov vectors are deflated
%
% -----------------------------------------------------------------------------
%
%  BANDITS: a Matlab Package of Band Krylov Subspace Iterations
%
%  Copyright (c) 2018-2019 Roland W. Freund
%  See LICENSE.txt for license
%
% -----------------------------------------------------------------------------

itmp = I.ph(1);
nv = norm(Vh(:,itmp));
%       
if n <= mc,
   if (d_f == 1) | (d_f == 2),
      def_fac = norm(R(:,n+m-mc),2);
   else
      def_fac = 1;
   end
else
   if (d_f == 1) | (d_f == 3),
      def_fac = nA;
   else
      def_fac = 1;   
   end
end
%        
act_defl_tol = def_fac * d_t;
%
if nv > act_defl_tol
   fv = 1;
else
   fprintf(' \n')
   disp('**--------------------------------------**')
%
   switch vw_f
      case 1
         disp('** Deflation of next candidate v vector **')
      case -1
         disp('** Deflation of next candidate w vector **')
      case 0
         disp('** Deflation of next candidate vector **')      
   end
%  
   disp('**--------------------------------------**')
   disp(['  Iteration index (n)    : ' num2str(n)])
   disp(['  Deflation Tolerance    : ' num2str(act_defl_tol)])
   disp(['  Norm of deflated vector: ' num2str(nv)])
%
   switch vw_f
      case 1   
         disp(['  New right block size   : ' num2str(mc-1)])
      case -1
         disp(['  New left block size    : ' num2str(mc-1)])
      case 0      
         disp(['  New block size         : ' num2str(mc-1)])    
     end
%  
   disp('**--------------------------------------**'), fprintf(' \n')
%
%  Shift the pointers of the candidate vectors
%
   for k = 1 : mc - 1,
      I.ph(k) = I.ph(k+1);
   end
%
%  Update I.nd, I.I, and I.pd
%
   if nv > 0,    
      I.nd = I.nd + 1;
      I.I = [I.I, n-mc];
      I.pd = [I.pd, itmp];
   else
      Vh(:,itmp) = zeros(size(Vh(:,itmp)));
   end
%
%  Update current block size and resize I.ph
%
   mc = mc - 1;
   I.ph = I.ph(1:mc);
%
end

