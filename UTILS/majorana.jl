#majorana utils
function M_OP(F :: F_OP)
	#transform Fermionic operator into Majorana operator
	Nmajs = 2*F.Nbods
	if Nmajs > 4
		error("Trying to build Majorana operator for more than 2-body fermionic operator, not implemented!")
	end
	if F.spin_orb
		error("Trying to build Majorana operator from spin-orb=true, not implemented!")
	end

	spin_orb = F.spin_orb
	filled_iso = zeros(Bool, Nmajs+1)
	filled_hetero = zeros(Bool, Nmajs+1)
	MBTS_ISO = Array{Float64}[]
	MBTS_HETERO = Array{Float64}[]
	body_sym = true
	N = F.N
	t_coeffs_iso = ones(Complex, Nmajs+1)
	t_coeffs_hetero = ones(Complex, Nmajs+1)

	#identity
	filled_iso[1] = true
	id_const = F.mbts[1][1]
	if F.filled[2]
		for i in 1:F.N
			id_const += F.mbts[2][i,i]
		end
	end
	if F.filled[3]
		for i in 1:F.N
			for j in 1:F.N
				id_const += F.mbts[3][i,i,j,j] + 0.5 * F.mbts[3][i,j,j,i]
			end
		end
	end
	push!(MBTS_ISO, [id_const])
	push!(MBTS_HETERO, [0])

	#one-body (i.e. 2-Majorana term)
	if F.Nbods > 0
		filled_iso[3] = true
		t_coeffs_iso[3] = 1im
		if F.filled[2]
			obt = copy(F.mbts[2])/2
		else
			obt = zeros(N, N)
		end

		if F.Nbods > 1 && F.filled[3]
			for i in 1:N
				for j in 1:N
					obt[i,j] += sum([F.mbts[3][i,j,k,k] for k in 1:N])
				end
			end
		end

		push!(MBTS_ISO, [0])
		push!(MBTS_ISO, obt)

		push!(MBTS_HETERO, [0])
		push!(MBTS_HETERO, [0])
	end

	#two-body
	if F.Nbods > 1
		push!(MBTS_ISO, [0])
		push!(MBTS_HETERO, [0])

		if F.filled[3] == false
			tbt_iso = [0]
			tbt_hetero = [0]
		else
			filled_iso[5] = true
			filled_hetero[5] = true
			tbt_hetero = -copy(F.mbts[3])/4
			tbt_iso = copy(F.mbts[3])

			for i in 1:N
				for k in 1:i-1
					for l in 1:N
						for j in 1:l-1
							tbt_iso[i,j,k,l] -= F.mbts[3][i,l,k,j]
						end
					end
				end
			end
			tbt_iso *= -0.5
		end
		
		push!(MBTS_ISO, tbt_iso)
		push!(MBTS_HETERO, tbt_hetero)
	end

	mbts_iso = tuple(MBTS_ISO...)
	mbts_hetero = tuple(MBTS_HETERO...)
	return M_OP(Nmajs, mbts_iso, mbts_hetero, t_coeffs_iso, t_coeffs_hetero, filled_iso, filled_hetero, body_sym, spin_orb, N)
end