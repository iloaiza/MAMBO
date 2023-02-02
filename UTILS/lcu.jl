#module for creating LCU from fragments and calculating the resulting 1-norm
function L1(F :: F_FRAG; debug = false, count = false)
	#calculates 1-norm of fragment
	#debug: run debugging routines
	if F.TECH == CSA()
		return CSA_L1(F, debug = debug, count = count)
	elseif F.TECH == DF()
		return DF_L1(F, debug = debug, count = count)
	elseif F.TECH == THC()
		return THC_L1(F, debug = debug, count = count)
	elseif F.TECH == OBF()
		return OBF_L1(F, debug = debug, count = count)
	else
		error("Trying to calculate LCU 1-norm decomposition for fermionic fragment with FRAGMENTATION_TECH=$(F.TECH), not implemented!")
	end
end

function SQRT_L1(F :: F_OP; count = false)
	#return minimal 1-norm for fermionic operator, does not scale well!
	of_OP = qubit_transform(to_OF(F))
	range = OF_qubit_op_range(of_OP)
	spec_range = (range[2] - range[1])/2

	if count
		return spec_range, 2
	else
		return spec_range
	end
end

function SQRT_L1(F :: F_FRAG; count = false)
	#lower bound of the 1-norm of a fragment, removes one-body correction
	Q_OP = qubit_transform(to_OF(F))
	obt_corr = ob_correction(F)
	Q_OP -= qubit_transform(to_OF(F_OP(obt_corr)))
	
	range = OF_qubit_op_range(Q_OP)
	spec_range = (range[2] - range[1])/2

	if count
		return spec_range, 2
	else
		return spec_range
	end
end

function CSA_L1(F :: F_FRAG; debug = false, count = false)
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

	if debug
		tbt = tbt_orb_to_so(cartan_2b_to_tbt(F.C))
		for i in 1:2F.N
			tbt[i,i,i,i] = 0
		end
		@show l1
		@show sum(abs.(tbt))/4
	end

	if count
		return l1, 2*(F.N^2)
	else
		return l1
	end
end

function DF_L1(F :: F_FRAG; debug = false, count = false)
	if F.spin_orb == false
		l1 = 0.5 * abs(F.coeff) * ((sum(abs.(F.C.λ)))^2)
		if count
			return l1, 1
		else
			return l1
		end
	else
		error("1-norm DF calculation not implemented for spin-orb=true")
	end
end

function THC_L1(F :: F_FRAG; debug = false, count = false)
	if F.spin_orb == false
		if count
			return abs.(F.coeff), 4
		else
			return abs(F.coeff)
		end
	else
		error("1-norm THC calculation not implemented for spin-orb=true")
	end
end

function OBF_L1(F :: F_FRAG; debug = false, count = false)
	if F.spin_orb == false
		l1 = sum(abs.(F.C.λ))
	else
		l1 = 0.5 * sum(abs.(F.C.λ))
	end

	if count
		return l1, F.N*(F.spin_orb+1)
	else
		return l1
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