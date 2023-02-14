# OPTIONS FOR WHAT ROUTINES SHOULD BE RAN
DO_CSA = true #perform Cartan Sub-Algebra (CSA) decomposition of Hamiltonian
DO_DF = true #perform Double-Factorization (DF) decomposition of Hamiltonian
DO_ΔE = true #obtain lower bound ΔE/2 of Hamiltonian, only for small systems!
DO_AC = true #do anticommuting grouping technique
DO_OO = true #do orbital optimization routine
DO_SQRT = true #obtain square-root factorization 1-norms
DO_MHC = true #do Majorana Hyper-Contraction routine
SYM_SHIFT = true #do symmetry-shift optimization routine
INT = true #do interaction picture optimization routines
verbose = false #verbose for sub-routines
COUNT = true #whether total number of unitaries should be counted
FBSS = true #whether fluid-body symmetry shift routine is done
DO_TROTTER = true #whether Trotter α is calculated

######## RUNNING CODE
mol_name = ARGS[1]

include("build.jl")

###### SAVELOAD ROUTINES FOR MOLECULAR HAMILTONIAN #######
FILENAME = DATAFOLDER*mol_name
if SAVING && isfile(FILENAME*".h5")
	fid = h5open(FILENAME*".h5", "cw")
	if haskey(fid, "MOLECULAR_DATA")
		println("Loading molecular data from $FILENAME.h5")
		MOL_DATA = fid["MOLECULAR_DATA"]
		h_const = read(MOL_DATA,"h_const")
		obt = read(MOL_DATA,"obt")
		tbt = read(MOL_DATA,"tbt")
		η = read(MOL_DATA,"eta")
		close(fid)
		H = F_OP((h_const,obt,tbt))
	else
		H, η = obtain_H(mol_name)
		println("""Saving molecular data in $FILENAME.h5 under group "MOLECULAR_DATA". """)
		if haskey(fid, "MOLECULAR_DATA")
			@warn "Trying to save molecular data to $FILENAME.h5, but MOLECULAR_DATA group already exists. Overwriting and migrating old file..."
			close(fid)
			oldfile(mol_name)
			fid = h5open(FILENAME*".h5", "cw")
		end
		create_group(fid, "MOLECULAR_DATA")
		MOL_DATA = fid["MOLECULAR_DATA"]
		MOL_DATA["h_const"] =  H.mbts[1]
		MOL_DATA["obt"] =  H.mbts[2]
		MOL_DATA["tbt"] =  H.mbts[3]
		MOL_DATA["eta"] =  η
		close(fid)
	end
else 
	H, η = obtain_H(mol_name)
	if SAVING
		println("""Saving molecular data in $FILENAME.h5 under group "MOLECULAR_DATA". """)
		fid = h5open(FILENAME*".h5", "cw")
		if haskey(fid, "MOLECULAR_DATA")
			@warn "Trying to save molecular data to $FILENAME.h5, but MOLECULAR_DATA group already exists."
			close(fid)
			oldfile(FILENAME)
			fid = h5open(FILENAME*".h5", "cw")
		end
		create_group(fid, "MOLECULAR_DATA")
		MOL_DATA = fid["MOLECULAR_DATA"]
		MOL_DATA["h_const"] =  H.mbts[1]
		MOL_DATA["obt"] =  H.mbts[2]
		MOL_DATA["tbt"] =  H.mbts[3]
		MOL_DATA["eta"] =  η
		close(fid)
	end
end

###### END: SAVELOAD ROUTINES FOR MOLECULAR HAMILTONIAN #######

RUN(H, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, DO_TROTTER=DO_TROTTER, DO_MHC = DO_MHC, COUNT = COUNT, verbose=verbose, name=FILENAME*".h5")

if SYM_SHIFT
	println("\n\nStarting symmetry-shift routine...")
	@time H_SYM, shifts = symmetry_treatment(H, verbose=verbose, SAVENAME=FILENAME*"_SYM.h5") # H = H_SYM + shifts[1]*Ne2 + shifts[2]*Ne
	println("Finished obtaining symmetry shifts, running routines for shifted Hamiltonian...")
	RUN(H_SYM, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, DO_TROTTER=DO_TROTTER, DO_MHC = DO_MHC, COUNT = COUNT, verbose=verbose, name=FILENAME*"_SYM.h5")
end

if INT
	println("\n\nStarting interaction picture routine...")
	@time H_INT = INTERACTION(H, SAVENAME=FILENAME*"_INT.h5")
	println("Finished obtaining interaction picture Hamiltonian, starting post-processing...")
	RUN(H_INT, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, DO_TROTTER=DO_TROTTER, DO_MHC = DO_MHC, COUNT = COUNT, verbose=verbose, name=FILENAME*"_INT.h5")
end

if FBSS
	println("\n\n Starting fluid-body symmetry shift routine...")
	println("FBSS optimization...")
	H_fbss = total_gsym_optimizer(H, η, verbose=verbose, SAVENAME=FILENAME*"_FBSS.h5")
	println("Running 1-norm routines...")
	RUN(H_fbss, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, DO_TROTTER=DO_TROTTER, DO_MHC = DO_MHC, COUNT = COUNT, verbose=verbose, name=FILENAME*"_FBSS.h5")
end

if FBSS && INT
	println("\n\n Starting interaction picture + fluid-body symmetry shift routines...")
	println("\nRunning before routine (H -> fbss -> int)")
	@time H_before = INTERACTION(H_fbss, SAVENAME=FILENAME*"_FBSS_INT.h5")
	RUN(H_before, DO_CSA = DO_CSA, DO_DF = DO_DF, DO_ΔE = DO_ΔE, DO_AC = DO_AC, DO_OO = DO_OO, DO_SQRT = DO_SQRT, DO_TROTTER=DO_TROTTER, DO_MHC = DO_MHC, COUNT = COUNT, verbose=verbose, name=FILENAME*"_FBSS_INT.h5")
end
