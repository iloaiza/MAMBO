mol_name = ARGS[1]

DO_CSA = true #perform Cartan Sub-Algebra (CSA) decomposition of Hamiltonian
DO_DF = true #perform Double-Factorization (DF) decomposition of Hamiltonian
DO_ΔE = true #obtain lower bound ΔE/2 of Hamiltonian, only for small systems!
DO_AC = true #do anticommuting grouping technique
DO_OO = true #do orbital optimization routine
DO_SQRT = true
SYM_SHIFT = true #do symmetry-shift optimization routine
INT = true #do interaction-picture optimization routines
verbose = false

@everywhere include("build.jl")
H, _ = obtain_H(mol_name)

RUN(H, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, verbose=verbose)

if SYM_SHIFT
	println("\n\nStarting symmetry-shift routine...")
	@time H_SYM, shifts = symmetry_treatment(H, verbose=verbose) # H = H_SYM + shifts[1]*Ne2 + shifts[2]*Ne
	println("Finished obtaining symmetry shifts, running routines for shifted Hamiltonian...")
	RUN(H_SYM, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, verbose=verbose)
end

if INT
	println("\n\nStarting interaction picture routine...")
	@time H_INT = INTERACTION(H)
	println("Finished obtaining interaction picture Hamiltonian, starting post-processing...")
	RUN(H_INT, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, verbose=verbose)
end