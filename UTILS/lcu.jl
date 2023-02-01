#module for creating LCU from fragments and calculating the resulting 1-norm
function L1(F :: F_FRAG; debug = false)
	#calculates 1-norm of fragment
	#debug: run debugging routines
	if F.TECH == CSA()
		return CSA_L1(F, debug = debug)
	elseif F.TECH == DF()
		return DF_L1(F, debug = debug)
	elseif F.TECH == THC()
		return THC_L1(F, debug = debug)
	elseif F.TECH == OBF()
		return OBF_L1(F, debug = debug)
	else
		error("Trying to calculate LCU 1-norm decomposition for fermionic fragment with FRAGMENTATION_TECH=$(F.TECH), not implemented!")
	end
end

function SQRT_L1(F :: F_OP)
	#return minimal 1-norm for fermionic operator, does not scale well!
	of_OP = qubit_transform(to_OF(F))
	range = OF_qubit_op_range(of_OP)
	return (range[2] - range[1])/2
end

function SQRT_L1(F :: F_FRAG; debug = false)
	#lower bound of the 1-norm of a fragment
	#debug: run debugging routines
	if F.TECH == CSA()
		return CSA_SQRT_L1(F, debug = debug)
	elseif F.TECH == DF()
		return DF_SQRT_L1(F, debug = debug)
	elseif F.TECH == THC()
		#THC already gives minimum 1-norm for a fragment since it only implements a rotated Z1
		return THC_L1(F, debug = debug)
	elseif F.TECH == OBF()
		#one-body fragments already give minimum 1-norm
		return OBF_L1(F, debug = debug)
	else
		error("Trying to calculate LCU square-root 1-norm decomposition for fermionic fragment with FRAGMENTATION_TECH=$(F.TECH), not implemented!")
	end
end

function CSA_L1(F :: F_FRAG; debug = false)
	if F.spin_orb == false
		idx = 0
		l1 = 0.0
		for i in 1:F.C.N
			for j in 1:i #sum i ≥ j
				idx += 1
				#λ_ij = λ_ji = F.C.λ[idx]
				if i != j
					#|λ_ij| + |λ_ji|
					l1 += 2*abs(F.C.λ[idx])
				else
					#|λ_ii|/2
					l1 += 0.5*abs(F.C.λ[idx])
				end
			end
		end
	else
		error("1-norm CSA calculation not implemented for spin-orb=true")
	end

	if debug == false
		return l1
	else
		tbt = tbt_orb_to_so(cartan_2b_to_tbt(F.C))
		@show l1
		for i in 1:2F.N
			tbt[i,i,i,i] = 0
		end
		@show sum(abs.(tbt))/4
		
		return l1
	end
end

function CSA_SQRT_L1(F :: F_FRAG; debug = false)
	if F.spin_orb
		error("Trying to calculate square-root 1-norm for spin-orb=true, not implemented!")
	end

	E_VALS = zeros(2^F.N)
	num_λ = length(F.C.λ)

	for i in 0:2^(F.N)-1
		n_arr = digits(i, base=2, pad=F.N)
		idx = 0
		for i1 in 1:F.N
			for i2 in 1:i1
				idx += 1
				if n_arr[i1] == 1 && n_arr[i2] == 1
					if i1 == i2
						E_VALS[i+1] += F.C.λ[idx]
					else
						E_VALS[i+1] += 2*F.C.λ[idx]
					end
				end
			end
		end
	end

	E_VALS *= 4

	if debug == false
		return (maximum(E_VALS) - minimum(E_VALS))/2
	else
		Q_OP = qubit_transform(to_OF(F))
		obt_corr = ob_correction(cartan_2b_to_tbt(F.C), F.spin_orb)
		obt_corr = obt_rotation(one_body_unitary(F.U[1]), corr_dbg)
		Q_OP -= qubit_transform(to_OF(F_OP(obt_corr)))

		return (maximum(E_VALS) - minimum(E_VALS))/2, Q_OP
	end
end

function DF_L1(F :: F_FRAG; debug = false)
	if F.spin_orb == false
		return 0.5 * abs(F.coeff) * ((sum(abs.(F.C.λ)))^2)
	else
		error("1-norm DF calculation not implemented for spin-orb=true")
	end
end

function DF_SQRT_L1(F :: F_FRAG; debug = false)
	if F.spin_orb == false
		if debug == false
			return 0.5 * abs(F.coeff) * ((sum(abs.(F.C.λ)))^2)
		else
			ϵ_pos = 0.5 * sum(F.C.λ + abs.(F.C.λ))
			ϵ_neg = 0.5 * abs(sum(F.C.λ - abs.(F.C.λ)))
			L1 = 0.5 * abs(F.coeff) * ((maximum([ϵ_pos, ϵ_neg]))^2)
			Q_OP = qubit_transform(to_OF(F))
			obt_corr = ob_correction(F)
			Q_OP -= qubit_transform(to_OF(F_OP(obt_corr)))
			return L1, Q_OP
		end
	else
		error("1-norm DF calculation not implemented for spin-orb=true")
	end
end

function THC_L1(F :: F_FRAG; debug = false)
	if F.spin_orb == false
		return abs(F.coeff)
	else
		error("1-norm THC calculation not implemented for spin-orb=true")
	end
end

function OBF_L1(F :: F_FRAG; debug = false)
	if F.spin_orb == false
		return sum(abs.(F.C.λ))
	else
		return 0.5 * sum(abs.(F.C.λ))
	end
end

function one_body_L1(H :: F_OP)
	#get one-body 1-norm after correction from 2-body term
	if H.spin_orb
		error("spin_orb = true, not implemented!")
	end
	obf = to_OBF(H.mbts[2] + ob_correction(H))

	return OBF_L1(obf)
end

function PAULI_L1(Q :: Q_OP)
	l1 = 0.0

	for pw in Q.paulis
		l1 += abs(pw.coeff)
	end

	return l1
end

function PAULI_L1(F :: F_OP)
	l1 = 0.0
	if F.Nbods > 2
		error("Trying to calculate Pauli cost for fermionic operator with more than 2-body terms, not implemented!")
	end
	if F.spin_orb
		error("Trying to calculate 1-norm for spin-orb=true, not implemented!")
	end

	if F.filled[2] && F.filled[3]
		obt_mod = F.mbts[2] + ob_correction(F)
		λ1 = sum(abs.(obt_mod))
	elseif F.filled[2]
		obt_mod = F.mbts[2]
		λ1 = sum(abs.(obt_mod))
	elseif F.filled[3]
		obt_mod = ob_correction(F)
		λ1 = sum(abs.(obt_mod))
	else
		λ1 = 0
	end
	
	if F.filled[3]
		λ2 = 0.5 * sum(abs.(F.mbts[3]))
		for r in 1:F.N
			for p in r+1:F.N
				for q in 1:F.N
					for s in q+1:F.N
						λ2 += abs(F.mbts[3][p,q,r,s] - F.mbts[3][p,s,r,q])
					end
				end
			end
		end
	else
		λ2 = 0
	end

	return λ1+λ2
end