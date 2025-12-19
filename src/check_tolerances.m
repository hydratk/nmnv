function [d_t,d_f,nA_f,nA,b_t,nJ_f,nJ] = check_tolerances(tol,n0)
%
%  This function is an auxiliary routine that extracts from the structure "tol"
%  various tolerances, flags, and norm estimates used to check for deflation
%  and breakdown.  It also checks if the deflation and breakdown tolerances are
%  in the range
%
%      (2^(-52) = ) eps <= tolerance < 1      (eps = machine precision)
%
%  and displays warnings in case they are not.  This auxiliary routine is
%  called in the functions "band_Arnoldi", "band_Lanczos", "sym_band_Lanczos",
%  "Herm_band_Lanczos", "Jsym_band_Lanczos", and "JHerm_band_Lanczos".
%
% -----------------------------------------------------------------------------
%
%  BANDITS: a Matlab Package of Band Krylov Subspace Iterations
%
%  Copyright (c) 2018-2019 Roland W. Freund
%  See LICENSE.txt for license
%
% -----------------------------------------------------------------------------

d_t = tol.defl_tol;
%
if (d_t < eps) || (d_t >= 1),
   fprintf(' \n')
   disp('**-------------------------------------**')    
   disp('** tol.defl_tol should be in the range **')
   disp('**                                     **')	    
   disp('**       eps <= tol.defl_tol < 1       **')
   disp('**                                     **')
   disp('**-------------------------------------**')
   fprintf(' \n') 
end
%
d_f = tol.defl_flag;
nA_f = tol.normA_flag;
%
if (nA_f == 1) & (n0 == 0),
   nA = 0;
else
   nA = tol.normA;
end
%
if nargout > 4
   b_t = tol.brk_tol;
   %
   if (b_t < eps) || (b_t >= 1),
      fprintf(' \n')
      disp('**------------------------------------**')    
      disp('** tol.brk_tol should be in the range **')
      disp('**                              	  **')	    
      disp('**       eps <= tol.brk_tol < 1       **')
      disp('**                                    **')
      disp('**------------------------------------**')
      fprintf(' \n') 
   end
end
%
if nargout > 5
   nJ_f = tol.normJ_flag;
   %
   if (nJ_f == 1) & (n0 == 0),
      nJ = 0;
   else
      nJ = tol.normJ;
   end
end

