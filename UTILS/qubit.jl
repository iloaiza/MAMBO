# qubit utils
function Q_OP(F :: F_OP, transformation = F2Q_map, tol = PAULI_TOL)
	return Q_OP(M_OP(F), transformation)
end

function majorana_pair_to_pauli(i, j, σ, n_qubits, transformation = F2Q_map)
	#transforms γiσ0*γjσ1 into Pauli word
	#σ ∈ {0,1}, corresponds to σ ∈ {false, true}
	if transformation != "jw" && transformation != "jordan-wigner"
		error("γiσ0*γjσ1 to Pauli word transformation not defined for $transformation!")
	end

	α = σ - 1
	bin_vec = zeros(Bool, 2*n_qubits)
	phase = -1im
	if i < j
		bin_vec[2i+α] = 1
		for n in 2i+α+1:2j+α-1
			bin_vec[n+n_qubits] = 1
		end
		bin_vec[2j+α] = 1
	elseif j < i
		bin_vec[2j+α] = 1
		bin_vec[2j+α + n_qubits] = 1
		for n in 2j+α+1:2i+α-1
			bin_vec[n+n_qubits] = 1
		end
		bin_vec[2i+α] = 1
		bin_vec[2i+α+n_qubits] = 1
	else
		bin_vec[2i+α+n_qubits] = 1
		phase *= -1
	end

	pw = pauli_word(bin_vec, phase)
	return pw
end

function Q_OP(M :: M_OP, transformation = F2Q_map, tol = PAULI_TOL)
	if M.spin_orb
		error("Transformation into Qubit operator not implemented for spin-orb = true!")
	end
	N = 2*M.N #number of qubits

	tol2 = tol^2

	n_paulis = Int(N^2/2 + N^4/4 + (N^2)*((N/2-1)^2)/16) #upper-bound for total number of Pauli words
	pws = [pw_zero(N) for i in 1:n_paulis]

	#identity term
	id_coeff = M.mbts_iso[1][1]


	curr_pauli = 1
	complex_flag = false

	#one-body terms
	if M.Nmajs ≥ 2
		if M.filled_hetero[1]
			error("Hetero index for Majorana is filled for 0-Majorana term, not defined!")
		end
		if M.filled_hetero[2]
			error("Hetero index for Majorana is filled for 1-Majorana term, not defined!")
		end
		if M.filled_hetero[3]
			error("Hetero index for Majorana is filled for 2-Majorana term, not defined!")
		end
		if M.filled_iso[2]
			error("Iso index for Majorana is filled for 1-Majorana term, qubit transform only defined for even terms!")
		end

		if M.filled_iso[3]
			for i in 1:M.N
				for j in 1:M.N
					ck = M.t_coeffs_iso[3] * M.mbts_iso[3][i,j]
					if abs2(ck) > tol2
						for σ in 0:1
							pws[curr_pauli] = ck * majorana_pair_to_pauli(i, j, σ, N, transformation)
							curr_pauli += 1
						end
					end
				end
			end
		end
	end

	#two-body terms
	if M.Nmajs ≥ 4
		if M.filled_hetero[4]
			error("Hetero index for Majorana is filled for 3-Majorana term, not defined!")
		end

		if M.filled_hetero[5]
			for i in 1:M.N
				for j in 1:M.N
					for k in 1:M.N
						for l in 1:M.N
							ck = M.t_coeffs_hetero[5] * M.mbts_hetero[5][i,j,k,l]
							if abs2(ck) > tol2
								for α in 0:1
									β = mod(α-1, 2)
									pw1 = majorana_pair_to_pauli(i, j, α, N, transformation)
									pw2 = majorana_pair_to_pauli(k, l, β, N, transformation)
									pws[curr_pauli] = ck * pw1 * pw2
									curr_pauli += 1
								end
							end
						end
					end
				end
			end
		end

		if M.filled_iso[4]
			error("Iso index for Majorana is filled for 3-Majorana term, qubit transform only defined for even terms!")
		end

		if M.filled_iso[5]
			for i in 1:M.N
				for l in 1:M.N
					for k in 1:i-1
						for j in 1:l-1
							ck = M.t_coeffs_iso[5] * M.mbts_iso[5][i,j,k,l]
							if abs2(ck) > tol2
								for α in 0:1
									pw1 = majorana_pair_to_pauli(i, j, α, N, transformation)
									pw2 = majorana_pair_to_pauli(k, l, α, N, transformation)
									pws[curr_pauli] = ck * pw1 * pw2
									curr_pauli += 1
								end
							end
						end
					end
				end
			end
		end
	end

	n_paulis = curr_pauli-1

	return Q_OP(N, n_paulis, id_coeff, pws[1:n_paulis])
end

function AC_group(Q :: Q_OP; ret_ops = false, verbose=false)
	is_grouped = zeros(Bool, Q.n_paulis)
	group_arrs = Array{Int64,1}[]
	vals_arrs = Array{Complex,1}[]

	vals_ord = [Q.paulis[i].coeff for i in 1:Q.n_paulis]
	ind_perm = sortperm(abs.(vals_ord))[end:-1:1]
	vals_ord = vals_ord[ind_perm]

	if verbose
		println("Running sorting-insertion algorithm")
		@show sum(vals_ord)
	end

	for i in 1:Q.n_paulis
    	if is_grouped[i] == false
    		curr_group = [i]
    		curr_vals = [vals_ord[i]]
    		is_grouped[i] = true
    		for j in i+1:Q.n_paulis
    			if is_grouped[j] == false
	    			if pws_is_anticommuting(Q.paulis[ind_perm[i]],Q.paulis[ind_perm[j]]) == 1
	    				antic_w_group = true
	    				for k in curr_group[2:end]
	    					if pws_is_anticommuting(Q.paulis[ind_perm[k]],Q.paulis[ind_perm[j]]) == 0
		    					antic_w_group = false
		    					break
		    				end
	    				end

	    				if antic_w_group == true
		    				push!(curr_group,j)
		    				push!(curr_vals,vals_ord[j])
		    				is_grouped[j] = true
		    			end
	    			end
	    		end
	    	end
    		push!(group_arrs,curr_group)
    		push!(vals_arrs, curr_vals)
    	end
    end

    if prod(is_grouped) == 0
    	println("Error, not all terms are grouped after AC-SI algorithm!")
    	@show is_grouped
    end

    num_groups = length(group_arrs)
    group_L1 = zeros(num_groups)
    for i in 1:num_groups
        for val in vals_arrs[i]
            group_L1[i] += abs2(val)
        end
    end

    L1_norm = sum(sqrt.(group_L1))

    if ret_ops == false
    	return L1_norm, num_groups
    else
    	OPS = Q_OP[]
    	for i in 1:num_groups
    		pws_i = [Q.paulis[ind_perm[j]] for j in group_arrs[i]]
    		q_op = Q_OP(Q.N, length(pws_i), 0, pws_i)
    		push!(OPS, q_op)
    	end

    	return L1_norm, OPS
    end
end

function AC_group(F :: F_OP; ret_ops = false, verbose=false)
	return AC_group(Q_OP(F), ret_ops=ret_ops, verbose=verbose)
end


