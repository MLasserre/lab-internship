function [y] = initial_data(r, r_0, v_left, v_right)
	if r <= r_0 then
		y = sign(v_left)*sqrt(1-(1-v_left^2)*(1-2/r)/(1-2/r_0))
	else
		y = sign(v_right)*sqrt(1-(1-v_right^2)*(1-2/r)/(1-2/r_0))
	end
endfunction
