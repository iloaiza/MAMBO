function INTERACTION(H)
	H0 = CSA_SD_greedy_decomposition(H :: F_OP, 1, verbose=false)[1]
	HR = H - H0

	return HR
end

function ORBITAL_OPTIMIZATION(H; verbose=true)
	H_rot, _, _ = orbital_l1_optimizer(H, verbose=verbose, ret_op=true)
	return H_rot
end

function RUN(H; DO_CSA = true, DO_DF = true, DO_ΔE = true, DO_AC = true, DO_OO = true, DO_SQRT = false, max_frags = 100, verbose=true)
	if DO_ΔE
		println("Obtaining 1-norm lower bound")
		@time λ_min = SQRT_L1(H)
		@show λ_min
	end

	if DO_CSA
		println("Doing CSA")
		max_frags = 100
		@time CSA_FRAGS = CSA_greedy_decomposition(H, max_frags, verbose=verbose)
	end

	if DO_DF
		println("\n\nDoing DF")
		@time DF_FRAGS = DF_decomposition(H, verbose=verbose)
	end

	println("\n\nCalculating 1-norms...")
	println("1-body:")
	@time λ1 = one_body_L1(H)
	@show λ1

	if DO_CSA
		println("\nCSA:")
		@time λ2_CSA = sum(L1.(CSA_FRAGS))
		@show (λ2_CSA, λ1 + λ2_CSA)
		if DO_SQRT
			println("Square-root routine...")
			t00 = time()
			l1 = SharedArray(zeros(length(CSA_FRAGS)))
			@sync @distributed for i in 1:length(CSA_FRAGS)
				l1[i] = SQRT_L1(CSA_FRAGS[i])
			end
			λ2_CSA_SQRT = sum(l1)
			println("Finished after $(time() - t00) seconds...")
			#@time λ2_CSA_SQRT = sum(SQRT_L1.(CSA_FRAGS))
			@show (λ2_CSA_SQRT, λ1 + λ2_CSA_SQRT)
		end
	end

	if DO_DF
		println("\nDF:")
		@time λ2_DF = sum(L1.(DF_FRAGS))
		@show (λ2_DF, λ1 + λ2_DF)
		#=
		if DO_SQRT
			println("Square-root routine...")
			@time λ2_DF_SQRT = sum(SQRT_L1.(DF_FRAGS))
			@show (λ2_DF_SQRT, λ1 + λ2_DF_SQRT)
		end
		# =#
	end

	println("\nPauli:")
	@time λPauli = PAULI_L1(H)
	@show λPauli
	
	if DO_AC
		println("\nAnti-commuting:")
		@time λAC, _ = AC_group(H, ret_ops=false)
		@show λAC
	end

	if DO_OO
		println("\nOrbital-rotation routine:")
		@time H_rot = ORBITAL_OPTIMIZATION(H, verbose=verbose)
		λOO_Pauli = PAULI_L1(H_rot)
		@show λOO_Pauli
		if DO_AC
			λOO_AC, _ = AC_group(H_rot, ret_ops=false)
			@show λOO_AC
		end
	end
end