
include("scp/scp_gusto.jl")
# include("scp/scp_trajopt.jl")
# include("scp/scp_mao.jl")

function add_constraint_category!(SCPCCat, func, dimtype, ind_time, ind_other...)
	SCPCCat[Symbol(func)] = ConstraintCategory(func, dimtype, ind_time, ind_other)
end
