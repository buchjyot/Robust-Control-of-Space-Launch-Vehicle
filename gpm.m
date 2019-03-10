function [gm_results,pm_results] = gpm(lp,indx)
% compute gain and phase margins

% lp = n x 2 complex matrix
% rows of lp correspond to n values of frequency
% 1st column of lp contains the complex gain for the loop transfer function
% 2nd column of lp contains the real frequency
% indx = integer for the lp that is stored with the results for convenience

% the results are stored in arrays with 3 columns
% where each row corresponds to a marginn a single frequency
% the number of rows can be 0, 1, 2, ...

% tmp = size(lp) ;
% n = tmp(1,1) - 1 ;
n = size(lp,1) ;

not_first = 0 ;

pm_results = [ ] ;
gm_results = [ ] ;

for i = 1:n 
   g = lp(i,1) ;
   w = lp(i,2) ;
   absg = abs(g) ;
   if( absg ~= 0 ) ,
      gain = 20 * log10(abs(g)) ;
      im = imag(g) ;
      re = real(g) ;
      if( not_first == 1 ) ,
         if( gain * gain_old <= 0 ) ,
            wc = w_old - ( w - w_old )/   ...
               ( gain - gain_old ) * gain_old ;
            gi = g_old + ( g - g_old )/   ...
               ( w - w_old )*( wc - w_old )  ;
            pm = atan2( -imag(gi) , -real(gi) ) * 180 / pi ;
            pm_results = [ pm_results ; indx wc pm ] ;
         end ;

         if( im * im_old <= 0 && re < 0 ) ,
            wp = w_old - ( w - w_old )/   ...
               ( im - im_old ) * im_old ;
            gi = g_old + ( g - g_old )/   ...
               ( w - w_old )*( wp - w_old )  ;
            gm = 20 * log10(abs(gi)) ;
            gm_results = [ gm_results ; indx wp gm ] ;
        end ;
      
      end ;
   
      not_first = 1 ;
      g_old = g ;
      w_old = w ;
      gain_old = gain ;
      im_old = im ;
   end ;   
   
end ;
clear g

