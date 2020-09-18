function [y] = solution_sigma(t_0,t_f,n,r_0,v_l,v_r,r_0,M)

	tt = linspace(t_0,t_f,n);

	function [z] = equation_sigma(r,v_l,v_r,r_0,M)
		for()
		z = function_R(r,v_l,v_r,r_0,M)..
           -function_R(r_0,v_l,v_r,r_0,M)
           -tt[i];
		end
	endfunction

	for i = t_0:n
		deff('[y]=g(x)','y = function_R(x, -0.48, -0.64, 10, 1) - function_R(10, -0.48, -0.64, 10, 1) - 0.01')
		fsolve(10, equation_sigma(r,v_l,v_r,r_0,M))
	end

endfunction
