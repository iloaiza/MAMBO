#utilities for Trotter

function trotter_α(FRAGS :: Array{F_FRAG})
	#calculates commutator norm ∑_{nm} ||[Hn,Hm]||
	#must include 1-body fragment as well
	M = length(FRAGS)
	α = 0.0

	for i in 1:M
		op_i = qubit_transform(to_OF(FRAGS[i]))
		for j in 1:i-1
			op_j = qubit_transform(to_OF(FRAGS[j]))
			comm_range = OF_qubit_op_range(of.commutator(op_i,op_j))
			α += 2*maximum(abs.(comm_range))
		end
	end

	return α
end

function trotter_β(FRAGS :: Array{F_FRAG})
	M = length(FRAGS)
	β = 0.0

	Δs = zeros(M)
	for i in 1:M
		curr_frag = FRAGS[i]
		Δs[i] = 2*SQRT_L1(curr_frag, count=false)
	end

	for i in 1:M
		for j in 1:i-1
			β += Δs[i] * Δs[j]
		end
	end

	return β
end
		


