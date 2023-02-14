function INTERACTION(H; SAVENAME=DATADIR*"INTERACTION.h5",verbose=false)
	H0 = CSA_SD_greedy_decomposition(H :: F_OP, 1, verbose=verbose, SAVENAME=SAVENAME)[1]
	HR = H - H0

	return HR
end

function ORBITAL_OPTIMIZATION(H; verbose=true, SAVELOAD=SAVING, SAVENAME=DATADIR*"OO.h5", do_Givens = OO_GIVENS)
	θmin = false
	if SAVELOAD
		fid = h5open(SAVENAME, "cw")
		if "ORBITAL_OPTIMIZATION" in keys(fid)
			OO = fid["ORBITAL_OPTIMIZATION"]
			θmin = read(OO, "theta")
		else
			create_group(fid, "ORBITAL_OPTIMIZATION")
			OO = fid["ORBITAL_OPTIMIZATION"]
		end
	end

	if θmin == false
		H_rot, _, _, θmin = orbital_l1_optimizer(H, verbose=verbose, ret_op=true)
		if SAVELOAD
			OO["theta"] = θmin
		end
	else
		println("Found saved θmin under filename $SAVENAME for orbital optimization, loaded orbital rotation...")
		if do_Givens
			U = givens_real_orbital_rotation(H.N, θmin)
		else
			U = real_orbital_rotation(H.N, θmin)
		end
		H_rot = F_OP_rotation(U, H)
	end

	if SAVELOAD
		close(fid)
	end

	return H_rot
end

function RUN(H; DO_CSA = true, DO_DF = true, DO_ΔE = true, DO_AC = true, DO_OO = true,
			 DO_SQRT = false, max_frags = 100, verbose=true, COUNT=false, DO_TROTTER = false,
			 DO_MHC = true, name = SAVING, SAVELOAD = SAVING)
	# Obtain 1-norms for different LCU methods. COUNT=true also counts number of unitaries in decomposition
	# CSA: Cartan sub-algebra decomposition
	# DF: Double Factorization
	# ΔE: Exact lower bound from diagonalization of H
	# AC: Anticommuting grouping
	# OO: Orbital rotation technique
	# SQRT: obtain square-root lower bound for non-optimal factorization methods (i.e. CSA)
	# TROTTER: obtain α upper-bound for Trotter error
	# MHC:L Majorana Hyper-Contraction
	# name: default name for saving, false means no saving is done

	if name == true
		name = DATADIR*"RUN.h5"
	end

	if DO_ΔE
		println("Obtaining 1-norm lower bound")
		if SAVELOAD
			fid = h5open(name, "cw")
			if haskey(fid, "dE")
				println("Found saved dE for file $name")
				λ_min = fid["dE"]
			else
				@time λ_min = SQRT_L1(H)
				fid["dE"] = λ_min
			end
			close(fid)
		end
		@show λ_min
	end

	if DO_TROTTER
		ob_frag = to_OBF(H.mbts[2] + ob_correction(H))
	end

	if DO_CSA
		println("Doing CSA")
		max_frags = 100
		@time CSA_FRAGS = CSA_greedy_decomposition(H, max_frags, verbose=verbose, SAVENAME=name)
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

	if DO_MHC
		println("\nMHC:")
		@time λ2_MHC_iter = iterative_schmidt(H.mbts[3], count=COUNT, tol=1e-6)
		@show λ1 + λ2_MHC_iter
		@time λ2_MHC_split = split_schmidt(H.mbts[3], count=COUNT, tol=1e-6)
		@show λ1 + λ2_MHC_split
	end

	if DO_CSA
		println("\nCSA:")
		@time λ2_CSA = sum(L1.(CSA_FRAGS, count=COUNT))
		@show λ1 + λ2_CSA
		if DO_SQRT
			println("Square-root routine...")
			@time λ2_CSA_SQRT = sum(SQRT_L1.(CSA_FRAGS, count=COUNT))
			@show λ1 + λ2_CSA_SQRT
		end
		if DO_TROTTER
			println("Starting Trotter routine for CSA...")
			CSA_TROT = CSA_FRAGS
			push!(CSA_TROT,ob_frag)
			@time β_CSA = trotter_β(CSA_TROT)
			@show β_CSA
		end
	end

	if DO_DF
		println("\nDF:")
		@time λ2_DF = sum(L1.(DF_FRAGS, count=COUNT))
		@show λ1 + λ2_DF
		if DO_TROTTER
			println("Starting Trotter routine for DF...")
			DF_TROT = DF_FRAGS
			push!(DF_TROT,ob_frag)
			@time β_DF = trotter_β(DF_TROT)
			@show β_DF
		end
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
		@time H_rot = ORBITAL_OPTIMIZATION(H, verbose=verbose, SAVENAME=name)
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