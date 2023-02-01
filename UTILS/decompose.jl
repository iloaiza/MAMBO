function CSA_greedy_decomposition(H :: F_OP, α_max; decomp_tol = ϵ, verbose=true)
	F_rem = copy(H) #fermionic operator, tracks remainder after removing found greedy fragments
	F_rem.filled[1:2] .= false #only optimize 2-body tensor

	if H.spin_orb == true
		error("CSA decomposition should be done with Hamiltonian represented in orbitals, not spin-orbitals!")
	end

	cartan_L = cartan_2b_num_params(H.N)
	unitary_L = real_orbital_rotation_num_params(H.N)

	tot_L = cartan_L + unitary_L

	Farr = F_FRAG[]

	curr_cost = L2_partial_cost(F_rem)
	if verbose
		println("Initial L2 cost is $curr_cost")
	end
	α_curr = 0

	while curr_cost > decomp_tol && α_curr < α_max
		α_curr += 1
		#x = [λ..., θ...] for cartan(λ) and U(θ)
		if verbose
			@time sol = CSA_greedy_step(F_rem)
			println("Current L2 cost after $α_curr fragments is $(sol.minimum)")
		else
			sol = CSA_greedy_step(F_rem)
		end
		frag = CSA_x_to_F_FRAG(sol.minimizer, H.N, H.spin_orb, cartan_L)
		push!(Farr, frag)
		F_rem = F_rem - to_OP(frag)
		curr_cost = sol.minimum
	end

	if verbose
		println("Finished CSA decomposition, total number of fragments is $α_curr, remainder L2-norm is $curr_cost")
	end

	if curr_cost > decomp_tol
		println("Warning, CSA decomposition did not converge, remining L2-norm is $curr_cost")
	end

	return Farr
end

function CSA_greedy_step(F :: F_OP, do_svd = SVD_for_CSA)
	cartan_L = cartan_2b_num_params(F.N)
	unitary_L = real_orbital_rotation_num_params(F.N)

	x0 = zeros(cartan_L + unitary_L)
	if do_svd
		frag = tbt_svd_1st(F.mbts[3]) #build initial guess from largest SVD fragment
		#frag = tbt_svd_avg(F.mbts[3]) #combine SVD fragments with largest unitary for initial guess
		x0[1:cartan_L] = frag.C.λ
		x0[cartan_L+1:end] = frag.U[1].θs
	else
		x0[cartan_L+1:end] = 2π*rand(unitary_L)
	end

	function cost(x)
		Fx = CSA_x_to_F_FRAG(x, F.N, F.spin_orb, cartan_L)
		return L2_partial_cost(F, to_OP(Fx))
	end

	return optimize(cost, x0, BFGS())
end

function CSA_SD_greedy_decomposition(H :: F_OP, α_max; decomp_tol = ϵ, verbose=true)
	#same as CSA decomposition, but includes optimization of one-body term
	F_rem = copy(H) #fermionic operator, tracks remainder after removing found greedy fragments
	F_rem.filled[1] = false #only optimize 1-body and 2-body tensors

	if H.spin_orb == true
		error("CSA_SD decomposition should be done with Hamiltonian represented in orbitals, not spin-orbitals!")
	end

	cartan_L = cartan_2b_num_params(H.N)
	unitary_L = real_orbital_rotation_num_params(H.N)

	tot_L = cartan_L + unitary_L + H.N

	Farr = F_FRAG[]

	curr_cost = L2_partial_cost(F_rem)
	if verbose
		println("Initial L2 cost is $curr_cost")
	end
	α_curr = 0

	while curr_cost > decomp_tol && α_curr < α_max
		α_curr += 1
		#x = [λ..., θ...] for cartan(λ) and U(θ)
		if verbose
			@time sol = CSA_SD_greedy_step(F_rem)
			println("Current L2 cost after $α_curr fragments is $(sol.minimum)")
		else
			sol = CSA_SD_greedy_step(F_rem)
		end
		frag = CSA_SD_x_to_F_FRAG(sol.minimizer, H.N, H.spin_orb, cartan_L)
		push!(Farr, frag)
		F_rem = F_rem - to_OP(frag)
		curr_cost = sol.minimum
	end

	if verbose
		println("Finished CSA_SD decomposition, total number of fragments is $α_curr, remainder L2-norm is $curr_cost")
		if curr_cost > decomp_tol
			println("Warning, CSA_SD decomposition did not converge, remining L2-norm is $curr_cost")
		end
	end

	return Farr
end

function CSA_SD_greedy_step(F :: F_OP, do_svd = SVD_for_CSA_SD)
	cartan_L = cartan_2b_num_params(F.N)
	unitary_L = real_orbital_rotation_num_params(F.N)

	x0 = zeros(cartan_L + unitary_L + F.N)
	if do_svd
		frag = tbt_svd_1st(F.mbts[3]) #build initial guess from largest SVD fragment
		#frag = tbt_svd_avg(F.mbts[3]) #combine SVD fragments with largest unitary for initial guess
		x0[F.N+1:F.N+cartan_L] = frag.C.λ
		x0[F.N+cartan_L+1:end] = frag.U[1].θs
	else
		x0[F.N+cartan_L+1:end] = 2π*rand(unitary_L)
	end

	function cost(x)
		Fx = CSA_SD_x_to_F_FRAG(x, F.N, F.spin_orb, cartan_L)
		return L2_partial_cost(F, to_OP(Fx))
	end

	return optimize(cost, x0, BFGS())
end

function THC_fixed_decomposition(Ftarg :: F_OP, α, θ0 = 2π*rand(Ftarg.N-1, α), ζ0 = zeros(Int(α*(α+1)/2)))
	#do THC decomp, will find angles corresponding to α rotations
	if size(θ0)[2] < α #starting θ guess does not cover all
		diff = α - size(θ0)[2]
		θ0 = hcat(θ0, 2π*rand(Ftarg.N-1, diff))
	end

	x0 = cat(ζ0, reshape(θ0, :), dims=1)

	function cost(x)
		FRAGSx = THC_x_to_F_FRAGS(x, α, Ftarg.N)
		F = to_OP(FRAGSx[1])
		for i in 2:length(FRAGSx)
			F += to_OP(FRAGSx[i])
		end

		return L2_partial_cost(Ftarg, F)
	end

	return optimize(cost, x0, BFGS())
end

function THC_iterative_decomposition(H :: F_OP, α_max, decomp_tol = ϵ)
	F = copy(H)
	F.filled[1:2] .= false

	println("Initial L2 cost is $(L2_partial_cost(F))")
	@time sol = THC_fixed_decomposition(F, 1)

	i = 2
	curr_cost = sol.minimum
	println("L2 cost with THC dimension 1 is $curr_cost")

	while i <= α_max && curr_cost > decomp_tol
		old_ζL = Int(i*(i-1)/2)
		ζ0 = zeros(Int(i*(i+1)/2))
		ζ0[1:Int(i*(i-1)/2)] = sol.minimizer[1:old_ζL]
		θ0 = reshape(sol.minimizer[old_ζL+1:end],H.N-1,i-1)
		@time sol = THC_fixed_decomposition(F, i, θ0, ζ0)
		i += 1
		curr_cost = sol.minimum

		println("L2 cost with THC dimension $(i-1) is $curr_cost")
	end

	return THC_x_to_F_FRAGS(sol.minimizer, i-1, H.N)
end


function DF_decomposition(H :: F_OP; tol=SVD_tol, tiny=SVD_tiny, verbose=true, debug=true, do_Givens=DF_GIVENS)
	#do double-factorization
	#do_Givens will try to transform each orbital rotation into Givens rotations, false returns one-body rotation matrices directly
	if verbose
		println("Starting Double-Factorization routine")
	end
	if H.spin_orb
		println("Doing Double-Factorization for spin-orb=true, be wary of results...")
	end
	
	n = H.N
	N = n^2


	tbt_full = reshape(H.mbts[3], (N,N))
	tbt_res = Symmetric(tbt_full)
	if sum(abs.(tbt_full - tbt_res)) > tiny
		println("Non-symmetric two-body tensor as input for SVD routine, calculations might have errors...")
		tbt_res = tbt_full
	end

	if verbose
		println("Diagonalizing two-body tensor")
		@time Λ,U = eigen(tbt_res)
	else
		Λ,U = eigen(tbt_res)
	end
	## U*Diagonal(Λ)*U' == tbt_res
	ind=sortperm(abs.(Λ))[end:-1:1]
    Λ = Λ[ind]
    U=U[:,ind]
	
	num_ops = N
    for i in 1:N
    	if abs(Λ[i]) < tol
    		if verbose
    			println("Truncating SVD for coefficients with magnitude smaller or equal to $(abs(Λ[i])), using $(i-1) fragments")
    		end
    		num_ops = i-1
    		break
    	end
	end
	Λ = Λ[1:num_ops]
	U = U[:,1:num_ops]

	FRAGS = [initialize_FRAG(n, DF()) for i in 1:num_ops]
	TBT = zeros(n,n,n,n)

	for i in 1:num_ops
        full_l = reshape(U[:, i], (n,n))
        cur_l = Symmetric(full_l)
        sym_dif = sum(abs.(cur_l - full_l))
        if sym_dif > tiny
        	if sum(abs.(full_l + full_l')) > tiny
				error("SVD operator $i is neither Hermitian or anti-Hermitian, cannot do double factorization into Hermitian fragment!")
			end
        	cur_l = Hermitian(1im * full_l)
        	Λ[i] *= -1
        end
    
		ωl, Ul = eigen(cur_l)
		if do_Givens
			if sum(abs.(imag.(log(Ul)))) > 1e-8
				Rl = f_matrix_rotation(n, Ul)
				C = cartan_1b(false, ωl, n)
				FRAGS[i] = F_FRAG(1, tuple(Rl), DF(), C, n, false, Λ[i], true)
			else
				Rl = SOn_to_MAMBO_full(Ul, verbose = false)
				if sum(abs2.(Ul - one_body_unitary(Rl))) > ϵ_Givens
					Rl = f_matrix_rotation(n, Ul)
				end
				C = cartan_1b(false, ωl, n)
				FRAGS[i] = F_FRAG(1, tuple(Rl), DF(), C, n, false, Λ[i], true)
			end
		else
			Rl = f_matrix_rotation(n, Ul)
			C = cartan_1b(false, ωl, n)
			FRAGS[i] = F_FRAG(1, tuple(Rl), DF(), C, n, false, Λ[i], true)
		end
	end

	if debug
		tbt_tot = zeros(n,n,n,n)
		for i in 1:num_ops
			tbt_tot += to_OP(FRAGS[i]).mbts[3]
		end

		L2_rem = sum(abs2.(tbt_tot - H.mbts[3]))
		if verbose
			@show L2_rem
		end
	end
    
    return FRAGS
end