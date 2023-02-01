# utilities for symplectic (binary) representation of Pauli words
struct pauli_bit
    bin :: Tuple{Bool,Bool}
end

function pauli_bit(i :: Int64)
    if i == 0
        return pauli_bit((0,0))
    elseif i == 1
        return pauli_bit((1,0))
    elseif i == 2
        return pauli_bit((1,1))
    elseif i == 3
        return pauli_bit((0,1))
    end
end

mutable struct pauli_word
    bits :: Array{pauli_bit,1} #array of single-qubit pauli matrices
    size :: Int64 #number of qubits in pauli word
    coeff :: Number #coefficient multiplying pauli word
end

function pauli_word(bits :: Array{pauli_bit}, coeff=1.0)
    pauli_size = length(bits)

    return pauli_word(bits,length(bits),coeff)
end

function pauli_num_from_bit(p :: pauli_bit)
    if p.bin[1] == false && p.bin[2] == false
        return 0
    elseif p.bin[1] == true && p.bin[2] == false
        return 1
    elseif p.bin[1] == true && p.bin[2] == true
        return 2
    elseif p.bin[1] == false && p.bin[2] == true
        return 3
    end
end

function ϵ_pauli(i,j)
    #return what k corresponds for Pauli multiplication and sign
    if i == j
        return 1, 0
    elseif i == 0
        return 1, j
    elseif j == 0
        return 1, i
    else
        if i == 1 && j == 2
            return 1im, 3
        elseif i == 2 && j == 1
            return -1im, 3
        elseif i == 2 && j == 3
            return 1im, 1
        elseif i == 3 && j == 2
            return -1im, 1
        elseif i == 3 && j == 1
            return 1im, 2
        elseif i == 1 && j == 3
            return -1im, 2
        end
    end
end


function bit_sum(b1 :: pauli_bit, b2 :: pauli_bit)
    #σi*σj
    i = pauli_num_from_bit(b1)
    j = pauli_num_from_bit(b2)
    phase, k = ϵ_pauli(i, j)

    return pauli_bit(k), phase
end

import Base.*

function *(a :: pauli_word, b :: pauli_word)
    len1 = a.size
    len2 = b.size

    if len1 != len2
        error("Unequal lengths!")
    end

    bits = pauli_bit[]
    sizehint!(bits, len1)

    tot_phase = 1 + 0im
    for i in 1:len1 
        prod_result, phase = bit_sum(a.bits[i], b.bits[i])
        push!(bits, prod_result)
        tot_phase *= phase
    end

    coeff = a.coeff * b.coeff * tot_phase

    return pauli_word(bits, coeff)
end

function *(pw :: pauli_word, b :: Number)
	return pauli_word(pw.bits, pw.size, pw.coeff * b)
end

function *(b :: Number, pw :: pauli_word)
	return pauli_word(pw.bits, pw.size, pw.coeff * b)
end

import Base.one

one(T::Type{pauli_bit}) = pauli_bit(0)

function pw_identity(n_qubits)
	bin = zeros(Bool, 2*n_qubits)
	return pauli_word(bin)
end

function bin_vec_to_pw(b_vec, n_qubits = Int(size(b_vec)[1]/2), coeff=1)
    bits = ones(pauli_bit, n_qubits)

    for i in 1:n_qubits
        bits[i] = pauli_bit((b_vec[i], b_vec[n_qubits+i]))
    end

    return pauli_word(bits, coeff)
end

function pauli_word(bin_vec :: Array{Bool}, coeff = 1)
	return bin_vec_to_pw(bin_vec, Int(size(bin_vec)[1]/2), coeff)
end

function pw_zero(n_qubits)
	bits = ones(pauli_bit, n_qubits)

	return pauli_word(bits, 0)
end

function pw_to_bin_vec(pw)
    n_qubits = pw.size
    b_vec = zeros(Bool, 2*n_qubits)

    for (i,bit) in enumerate(pw.bits)
        b_vec[i] = bit.bin[1]
        b_vec[i+n_qubits] = bit.bin[2]
    end

    return b_vec
end

function bin_vec_prod(bin1,bin2,n_qubits=Int(length(bin1)/2))
    #returns, up to a phase, pauli corresponding to bin1*bin2
    bin3 = zeros(Bool,2*n_qubits)
    bin3 .= (bin1 + bin2) .% 2
end

function binary_is_anticommuting(bin1, bin2, n_qubits)
	#check if two binary vectors (i.e. Pauli words) are anticommuting
    return Bool(sum(bin1[1:n_qubits] .* bin2[n_qubits+1:end] + bin1[n_qubits+1:end] .* bin2[1:n_qubits]) % 2)
end

function pws_is_anticommuting(pw1, pw2)
	bin1 = pw_to_bin_vec(pw1)
	bin2 = pw_to_bin_vec(pw2)

	return binary_is_anticommuting(bin1, bin2, pw1.size)
end
