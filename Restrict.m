function [A_list, R_list,max_level] =Restrict(A, direct_N)
N =size(A,1);
n =round(sqrt(N));

A_list={};
R_list={};

level=1;
A_list(level) = {A};

while N> direct_N

    % 		coarse_n= floor((n - 1) / 2);
    % 		coarse_N= coarse_n * coarse_n;
    %
    % 		R = sparse(coarse_N, N);
    % 		k = 0;
    % 		for j = 2 : 2 : n
    % 			for i = 2 : 2 : n
    % 				k = k + 1;
    % 				fine_grid_k=(j - 1) * n + i;
    %
    % 				R(k,fine_grid_k - n- 1) = 0.0625;
    % 				R(k,fine_grid_k - n )  = 0.125;
    % 				R(k,fine_grid_k - n+1) = 0.0625;
    %
    % 				R(k,fine_grid_k - 1)   = 0.125;
    % 				R(k,fine_grid_k    )   = 0.25;
    % 				R(k,fine_grid_k + 1)   = 0.125;
    %
    % 				R(k,fine_grid_k + n - 1) = 0.0625;
    % 				R(k,fine_grid_k + n) = 0.125;
    % 				R(k,fine_grid_k + n + 1) = 0.0625;
    % 			end
    % 		end
    %
    % 		R_list(level) = {R};
    %
    % 		P = 4 * R';
    % 		A = R * A * P;
    %
    % 		level = level + 1;
    % 		A_list(level) = {A};
    % 		n = coarse_n;
    % 		N = coarse_N;
    % 	end
    % 	max_level = level;
    %
    % end




    nc = (n - 1)/2;          % coarse grid size
    Nc = nc^2;

    R = sparse(Nc, N);

    kc = 0;
    for jc = 1:nc
        for ic = 1:nc
            kc = kc + 1;

            i = 2*ic;
            j = 2*jc;
            kf = (j-1)*n + i;

            % 3x3 full-weighting stencil
            R(kc, kf-n-1) = 1/16;
            R(kc, kf-n  ) = 1/8;
            R(kc, kf-n+1) = 1/16;

            R(kc, kf-1  ) = 1/8;
            R(kc, kf    ) = 1/4;
            R(kc, kf+1  ) = 1/8;

            R(kc, kf+n-1) = 1/16;
            R(kc, kf+n  ) = 1/8;
            R(kc, kf+n+1) = 1/16;
        end
    end

    R_list{level} = R;
    P = 4*R';                 % linear interpolation
    A = R*A*P;            
    level = level + 1;
    A_list{level} = A;

    n = nc;
    N = Nc;
end

max_level = level;
end