#routines for finding fluid-body symmetry shifts
function gsym_params_to_F_OP(λvec, θvec, num_elecs, s=1)
	#builds two-body tensor corresponding to s2 = (O*Ne + Ne*O)/2
	#O = U'(θ) one-body-pol(λ) U(θ)
	#also returns s1 = η*t2*O
	#note t2 can be set = 1 and full flexibility remains in λvec
	U = givens_real_orbital_rotation(length(λvec), θvec)
	C = cartan_1b(false, λvec)

	#U_mat = one_body_unitary(U)
	obt = s*cartan_obt_rotation(U,cartan_1b_to_obt(C))
	tbt = zeros(U.N,U.N,U.N,U.N)

	for i in 1:U.N
		for j in 1:U.N
			for k in 1:U.N
				tbt[i,j,k,k] += obt[i,j]
				tbt[i,i,j,k] += obt[j,k]
			end
		end
	end

	tbt *= 0.5

	return F_OP(([0],num_elecs*obt,tbt))
end

function gsym_optimizer(F :: F_OP, num_elecs; s=true, verbose=false)
	#find best parameters for generalized coefficient shift
	#include t2 for also optimizing obt*Ne
	#does not include Ne and Ne^2 optimization
	if F.spin_orb
		error("Generalized symmetry shift not defined for spin-orb=true!")
	end

	nθ = real_orbital_rotation_num_params(F.N)
	x0 = zeros(F.N + nθ + s)

	function cost(x, s_val = s)
		if s_val
			Sx = gsym_params_to_F_OP(x[2:F.N+1], x[F.N+2:end], num_elecs, x[1])
		else	
			Sx = gsym_params_to_F_OP(x[1:F.N], x[F.N+1:end], num_elecs)
		end

		return PAULI_L1(F-Sx)
	end

	if verbose
		println("Starting 1-norm cost:")
		@show cost(x0)
		@time sol = optimize(cost, x0, BFGS())
		println("Final 1-norm cost:")
		@show sol.minimum
	else
		sol = optimize(cost, x0, BFGS())
	end

	if s
		S = gsym_params_to_F_OP(sol.minimizer[2:F.N+1], sol.minimizer[F.N+2:end], num_elecs, sol.minimizer[1])
	else
		S = gsym_params_to_F_OP(sol.minimizer[1:F.N], sol.minimizer[F.N+1:end], num_elecs)
	end

	return F - S
end

function total_gsym_optimizer(F :: F_OP, num_elecs; s=true, verbose=false, SAVELOAD = SAVING, SAVENAME=DATAFOLDER*"FBSS.h5")
	#find best parameters for generalized coefficient shift
	#includes optimization with Ne and Ne^2
	#s=true includes a coefficient multypling obt. This coefficient is redundant since it appears in all λs, but sometimes gives better convergence
	if F.spin_orb
		error("Generalized symmetry shift not defined for spin-orb=true!")
	end

	nθ = real_orbital_rotation_num_params(F.N)
	x0 = zeros(F.N + nθ + 2 + s)
	if s
		x0[3] = 1
	end

	Ne, Ne2 = symmetry_builder(F)


	x = false

	if SAVELOAD
		fid = h5open(SAVENAME, "cw")
		if haskey(fid, "FBSS")
			FBSS_group = fid["FBSS"]
			if haskey(FBSS_group, "shifts")
				x = FBSS_group["shifts"]
				close(fid)
			end
		else
			create_group(fid, "FBSS")
			FBSS_group = fid["FBSS"]
		end
	end

	if SAVELOAD == false || x == false
		function cost(x, s_val = s)
			if s_val
				Sx = gsym_params_to_F_OP(x[4:F.N+3], x[F.N+4:end], num_elecs, x[3])
			else
				Sx = gsym_params_to_F_OP(x[3:F.N+2], x[F.N+3:end], num_elecs)
			end
			Sx += x[1]*Ne + x[2]*Ne2
			return PAULI_L1(F-Sx)
		end

		if verbose
			println("Starting 1-norm cost:")
			@show cost(x0)
			@time sol = optimize(cost, x0, BFGS())
			println("Final 1-norm cost:")
			@show sol.minimum
		else
			sol = optimize(cost, x0, BFGS())
		end

		if SAVELOAD
			FBSS_group["shifts"] = sol.minimizer
			close(fid)
		end

		if s
			S = gsym_params_to_F_OP(sol.minimizer[4:F.N+3], sol.minimizer[F.N+4:end], num_elecs, sol.minimizer[3])
		else
			S = gsym_params_to_F_OP(sol.minimizer[3:F.N+2], sol.minimizer[F.N+3:end], num_elecs)
		end
		S += sol.minimizer[1]*Ne + sol.minimizer[2]*Ne2
	else
		println("Loaded FBSS parameters!")
		if s
			S = gsym_params_to_F_OP(x[4:F.N+3], x[F.N+4:end], num_elecs, x[3])
		else
			S = gsym_params_to_F_OP(x[3:F.N+2], x[F.N+3:end], num_elecs)
		end
		S += x[1]*Ne + x[2]*Ne2
	end

	return F - S
end

function naive_fbss(F :: F_OP, η)
	#return fbss routine where O*Ne term is obtained directly from 2-body tensor
	if F.spin_orb
		error("Not defined for spin-orb = true!")
	end

	Ne, Ne2 = symmetry_builder(F)
	x0 = [1]

	function cost(x)
		obt = F.mbts[2]
		tbt = F.mbts[3]
		obt_corr = 2*x[1]*sum([tbt[:,:,k,k] for k in 1:F.N])
		for k in 1:F.N
			tbt[:,:,k,k] -= x[1]*F.mbts[3][:,:,k,k]
			tbt[k,k,:,:] -= x[1]*F.mbts[3][:,:,k,k]
		end
		F_curr = F_OP(([0],obt - η*obt_corr, tbt))

		return PAULI_L1(F - F_curr)
	end

	@show cost(x0)
	sol = optimize(cost, x0, BFGS())
	@show sol.minimum

	obt = F.mbts[2]
	tbt = F.mbts[3]
	obt_corr = 2*sol.minimizer[1]*sum([tbt[:,:,k,k] for k in 1:F.N])
	for k in 1:F.N
		tbt[:,:,k,k] -= sol.minimizer[1]*F.mbts[3][:,:,k,k]
		tbt[k,k,:,:] -= sol.minimizer[1]*F.mbts[3][:,:,k,k]
	end

	H_curr = F_OP(([0],obt - η*obt_corr, tbt))
	H_SYM, _ = symmetry_treatment(H_curr)

	@show PAULI_L1(H_SYM)

	error("Function not working!")
	return H_SYM
end

function l2_total_gsym_optimizer(F :: F_OP, num_elecs)
	#find best parameters for generalized coefficient shift
	#includes optimization with Ne and Ne^2
	if F.spin_orb
		error("Generalized symmetry shift not defined for spin-orb=true!")
	end

	nθ = real_orbital_rotation_num_params(F.N)
	x0 = zeros(F.N + nθ + 2)

	Ne, Ne2 = symmetry_builder(F)

	function cost(x)
		Sx = gsym_params_to_F_OP(x[3:F.N+2], x[F.N+3:end], num_elecs)
		Sx += x[1]*Ne + x[2]*Ne2
		return L2_total_cost(F-Sx)
	end

	@show cost(x0)
	sol = optimize(cost, x0, BFGS())
	@show sol.minimum

	S = gsym_params_to_F_OP(sol.minimizer[3:F.N+2], sol.minimizer[F.N+3:end], num_elecs)
	S += sol.minimizer[1]*Ne + sol.minimizer[2]*Ne2

	return F - S
end

function fbss_linprog(F :: F_OP, η; model="ipopt", verbose=true)
	if F.spin_orb
		error("Generalized symmetry shift not defined for spin-orb=true!")
	end

    if model == "highs"
        L1_OPT = Model(HiGHS.Optimizer)
    elseif model == "ipopt"
        L1_OPT = Model(Ipopt.Optimizer)
    else
        error("Not defined for model = $model")
    end
    
    if verbose == false
        set_silent(L1_OPT)
    end

    ν1_len = F.N^2
    ν2_len = F.N^4
    ν3_len = Int((F.N*(F.N-1)/2)^2)
    
    @variables(L1_OPT, begin
        t[1:2]
        obt[1:ν1_len]
        tbt1[1:ν2_len]
        tbt2[1:ν3_len]
        omat[1:ν1_len]
    end)

    @objective(L1_OPT, Min, sum(obt)+sum(tbt1)+sum(tbt2))

    obt_corr = ob_correction(F)
    #1-body 1-norm
    λ1 = zeros(ν1_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx += 1
    		λ1[idx] = F.mbts[2][i,j] + obt_corr[i,j]
    	end
    end

    τ_11 = zeros(ν1_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx += 1
    		if i == j
    			τ_11[idx] = 2*F.N
    		end
    	end
    end
    τ_12 = zeros(ν1_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx += 1
    		if i == j
    			τ_11[idx] = 1
    		end
    	end
    end
    T1 = zeros(ν1_len,ν1_len)
    T1 += Diagonal((η - F.N/2)*ones(ν1_len))
    idx1 = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx1 += 1
    		idx2 = 0
    		for k in 1:F.N
    			for l in 1:F.N
    				idx2 += 1
    				if i == j && k == l
 	   					T1[idx1,idx2] -= 1
 	   				end
 	   			end
 	   		end
 	   	end
 	end
 	@constraint(L1_OPT, low_1, λ1 - τ_11*t[1] - τ_12*t[2] + T1*omat - obt .<= 0)
 	@constraint(L1_OPT, high_1, λ1 - τ_11*t[1] - τ_12*t[2] + T1*omat + obt .<= 0)

 	#2-body αβ/βα 1-norm
 	λ2 = zeros(ν2_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:F.N
    			for l in 1:F.N
    				idx += 1
    				λ2[idx] = 0.5 * F.mbts[3][i,j,k,l]
    			end
    		end
    	end
    end

    τ_21 = zeros(ν2_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:F.N
    			for l in 1:F.N
    				idx += 1
    				if i == j && k == l
    					τ_21[idx] = 0.5
    				end
    			end
    		end
    	end
    end

    T2 = zeros(ν2_len,ν1_len)
    idx = 0
    idx_ij = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx_ij += 1
    		idx_kl = 0
    		for k in 1:F.N
    			for l in 1:F.N
    				idx += 1
    				idx_kl += 1
    				if i == j
    					T2[idx,idx_kl] += 1
    				end
    				if k == l
    					T2[idx,idx_ij] += 1
    				end
    			end
    		end
    	end
    end
    @constraint(L1_OPT, low_2, λ2 - τ_21*t[1] - 0.5*T2*omat - tbt1 .<= 0)
    @constraint(L1_OPT, high_2, λ2 - τ_21*t[1] - 0.5*T2*omat + tbt1 .<= 0)

    #2-body αα/ββ 1-norm
    λ3 = zeros(ν3_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:i-1
    			for l in 1:j-1
    				idx += 1
    				λ3[idx] = F.mbts[3][i,j,k,l] - F.mbts[3][i,l,k,j]
    			end
    		end
    	end
    end
    
    τ_31 = zeros(ν3_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:i-1
    			for l in 1:j-1
    				idx += 1
    				if i == j && k == l
    					τ_31[idx] += 1
    				end
    				if i == l && k == j
    					τ_31[idx] -= 1
    				end
    			end
    		end
    	end
    end
    
    T_dict = zeros(Int64,F.N,F.N)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		idx += 1
    		T_dict[i,j] = idx
    	end
    end

    T3 = zeros(ν3_len,ν1_len)
    idx = 0
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:i-1
    			for l in 1:j-1
    				idx += 1
    				
    				idx_ij = T_dict[i,j]
    				if k == l
    					T3[idx,idx_ij] += 1
    				end
    				
    				idx_kl = T_dict[k,l]
    				if i == j
    					T3[idx,idx_kl] += 1
    				end

    				idx_il = T_dict[i,l]
    				if k == j
    					T3[idx,idx_il] -= 1
    				end

    				idx_kj = T_dict[k,j]
    				if i == l
    					T3[idx,idx_kj] -= 1
    				end
    			end
    		end
    	end
    end

    @constraint(L1_OPT, low_3, λ3 - τ_31*t[1] - 0.5*T3*omat - tbt2 .<= 0)
    @constraint(L1_OPT, high_3, λ3 - τ_31*t[1] - 0.5*T3*omat + tbt2 .<= 0)

    optimize!(L1_OPT)

    t_opt = value.(t)
    o_opt = value.(omat)
    O = zeros(F.N,F.N)
    for i in 1:F.N
    	for j in 1:F.N
    		O[i,j] = o_opt[T_dict[i,j]]
    	end
    end

    Ne,Ne2 = symmetry_builder(F)
    
    s2_tbt = t_opt[1] * Ne2.mbts[3]
    for i in 1:F.N
    	for j in 1:F.N
    		for k in 1:F.N
    			s2_tbt[i,j,k,k] += 0.5*O[i,j]
    			s2_tbt[k,k,i,j] += 0.5*O[i,j]
    		end
    	end
    end
    s2 = F_OP(([0],[0],s2_tbt))

    s1_obt = t_vec[2]*Ne.mbts[2] - η*O
    s1 = F_OP(([0],s1_obt))

    return F - s2 - s1, s1+s2
end


