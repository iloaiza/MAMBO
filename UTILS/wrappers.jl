function INTERACTION(H)
	H0 = CSA_SD_greedy_decomposition(H :: F_OP, 1, verbose=false)[1]
	HR = H - H0

	return HR
end

function ORBITAL_OPTIMIZATION(H; verbose=true)
	H_rot, _, _ = orbital_l1_optimizer(H, verbose=verbose, ret_op=true)
	return H_rot
end

function RUN(H; DO_CSA = true, DO_DF = true, DO_ΔE = true, DO_AC = true, DO_OO = true, DO_SQRT = false, max_frags = 100, verbose=true, COUNT=false)
	# Obtain 1-norms for different LCU methods. COUNT=true also counts number of unitaries in decomposition
	# CSA: Cartan sub-algebra decomposition
	# DF: Double Factorization
	# ΔE: Exact lower bound from diagonalization of H
	# AC: Anticommuting grouping
	# OO: Orbital rotation technique
	# SQRT: obtain square-root lower bound for non-optimal factorization methods (i.e. CSA)

	if DO_ΔE
		println("Obtaining 1-norm lower bound")
		@time λ_min = SQRT_L1(H)
		@show λ_min
	end

	if DO_CSA
		println("Doing CSA")
		max_frags = 100
		@time CSA_FRAGS = CSA_greedy_decomposition(H, max_frags, verbose=verbose)
		println("Finished CSA decomposition for 2-body term using $(length(CSA_FRAGS)) fragments")
	end

	if DO_DF
		println("\n\nDoing DF")
		@time DF_FRAGS = DF_decomposition(H, verbose=verbose)
		println("Finished DF decomposition for 2-body term using $(length(DF_FRAGS)) fragments")
	end

	println("\n\nCalculating 1-norms...")
	println("1-body:")
	@time λ1 = one_body_L1(H, count=COUNT)
	@show λ1

	if DO_CSA
		println("\nCSA:")
		@time λ2_CSA = sum(L1.(CSA_FRAGS, count=COUNT))
		@show λ1 + λ2_CSA
		if DO_SQRT
			println("Square-root routine...")
			@time λ2_CSA_SQRT = sum(SQRT_L1.(CSA_FRAGS, count=COUNT))
			@show λ1 + λ2_CSA_SQRT
		end
	end

	if DO_DF
		println("\nDF:")
		@time λ2_DF = sum(L1.(DF_FRAGS, count=COUNT))
		@show λ1 + λ2_DF
	end

	println("\nPauli:")
	@time λPauli = PAULI_L1(H, count=COUNT)
	@show λPauli
	
	if DO_AC
		println("\nAnti-commuting:")
		@time λAC, N_AC = AC_group(H, ret_ops=false)
		if COUNT
			@show λAC, N_AC
		else
			@show λAC
		end
	end

	if DO_OO
		println("\nOrbital-rotation routine:")
		@time H_rot = ORBITAL_OPTIMIZATION(H, verbose=verbose)
		λOO_Pauli = PAULI_L1(H_rot, count=COUNT)
		@show λOO_Pauli
		if DO_AC
			λOO_AC, N_OO_AC = AC_group(H_rot, ret_ops=false)
			if COUNT
				@show λOO_AC, N_OO_AC
			else
				@show λOO_AC
			end
		end
	end
end