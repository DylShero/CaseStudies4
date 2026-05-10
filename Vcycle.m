function x = Vcycle(level, A_list, R_list, b, x0, direct_n, weight, pre_smooth, pos_smooth)

% level      Current level (initial=1)
% A_list     Array of coefficient matrices on each level
% R_list     Array of restriction operators on each level
% b          Right-hand side
% x0         Initial Guess
% direct_n   Threshold for directly solving Ax = b
% weight     Scaling coefficient between restriction and prolongation
%            operators
% pre_steps  Number of iters of pre-smoothing
% pos_steps  Number of iters of post-smoothing


	
	% Load coefficient matrix
	A = cell2mat(A_list(level));
	
	%Solve directly if problem is small enough
	n =size(b, 1);
	if (n <= direct_n)
		x = A \ b;
		return;
    end

    x = x0

    w = 2/3;
    D = diag(A);

	% Pre-smoothing
	for i = 1:pre_smooth
        r = b - A * x;
        x = x + w * (r ./ D); 
    end
	
	% Load restriction operator and construct interpolation operator
	R = cell2mat(R_list(level));
	P =  R' * weight;       % makes sure the coarse grid has the correct scaling when interpolated back to the fine grid. 
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	r = b - A * x;
    r_H = R * r;
	
	% Implement the V-cycle recursively
	x0  = zeros(coarse_n, 1);
	e = Vcycle(level + 1, A_list, R_list, r_H, x0, direct_n, weight, pre_smooth, pos_smooth);
	
	% Transfer error to fine grid and correct
	x = x + P * e;
	
	% Post-smoothing
	for i = 1:pos_smooth
        r = b - A * x;
        x = x + w * (r ./ D); 
    end
     
end