function result = save_result(aflag,V,Vh,mc,Iv,T,rho,tol,D,W,Wh,pc,Iw,Tt,eta)
%
%  This function is an auxiliary routine that is used to save the values of
%  variables in the output structure "result".  It is called in the functions
%  "band_Arnoldi", "band_Lanczos", "sym_band_Lanczos", "Herm_band_Lanczos",
%  "Jsym_band_Lanczos", and "JHerm_band_Lanczos".
%
%  The input flag "aflag" is used to save the correct quantities for the
%  various algorithms:
%
%  aflag = 0   "band_Arnoldi"
%
%  aflag = 1   "band_Lanczos" 
%
%  aflag = 2   "Herm_band_Lanczos"
%
%  aflag = 3   "sym_band_Lanczos", "Jsym_band_Lanczos", "JHerm_band_Lanczos"
%
% -----------------------------------------------------------------------------
%
%  BANDITS: a Matlab Package of Band Krylov Subspace Iterations
%
%  Copyright (c) 2018-2019 Roland W. Freund
%  See LICENSE.txt for license
%
% -----------------------------------------------------------------------------

switch aflag
   case 0
      result.V = V;
      result.Vh_defl = Vh;
      result.mc = mc;
      result.Iv = Iv;
      result.rho = rho;
      result.H = T;
      result.Iv = Iv;
   case 1
      result.V = V;
      result.Vh_defl = Vh;
      result.mc = mc;
      result.Iv = Iv;
      result.rho = rho;
      result.T = T;
%      
      result.W = W;
      result.Wh_defl = Wh;
      result.pc = pc;
      result.Iw = Iw;
      result.eta = eta;
      result.Tt = Tt;
%
      result.D = D;
   case 2
      result.V = V;
      result.Vh_defl = Vh;
      result.mc = mc;
      result.Iv = Iv;
      result.rho = rho;
      result.T = T;  
   case 3
      result.V = V;
      result.Vh_defl = Vh;
      result.mc = mc;
      result.Iv = Iv;
      result.rho = rho;
      result.T = T;  
%    
      result.D = D;
end
%  
result.tol = tol;
      
