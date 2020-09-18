function [y] = function_R(r, v_l, v_r, r_0, M)

	K_l = (1-v_l^2)/(1-2*M/r_0)
	K_r = (1-v_r^2)/(1-2*M/r_0)

	y = 2*sign(v_l + v_r)/(K_r - K_l)..
       *(..
          (r-4*M)/(r-2*M)..
         *(..
            sqrt(1-K_l*(1-2*M/r))-sqrt(1-K_r*(1-2*M/r))..
          )..
         -M*(4-3*K_r)..
         *(..
            log(..
                 r*sqrt(1-K_l)*sqrt(1-K_r*(1-2*M/r))..
                +(M-r)*K_r + r..
               )..
           +log(..
                 2*r/(r-2*M)*sqrt(1-K_r*(1-2*M/r)) - K_r..
               )..
           )..
         +M*(4-3*K_l)..
         *(..
            log(..
                 r*sqrt(1-K_l)*sqrt(1-K_l*(1-2*M/r)) + (M-r)*K_l + r..
               )..
           -log(2*r/(r-2*M)*sqrt(1-K_l*(1-2*M/r)) - K_l..
               )..
          )..
         )

endfunction

