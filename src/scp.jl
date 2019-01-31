
include("scp/scp_gusto.jl")
# include("scp/scp_trajopt.jl")
# include("scp/scp_mao.jl")

function add_constraint_category!(SCPCCat, func, dimtype::Symbol, ind_time, ind_other...)
	SCPCCat[Symbol(func)] = ConstraintCategory(func, dimtype, ind_time, ind_other)
end

function add_constraint_category!(SCPCCat, func, goal::G, dimtype::Symbol, ind_other...) where G <: Goal{T} where T
	SCPCCat[Symbol(func)] = ConstraintCategory(func, dimtype, goal.ind_coordinates, ind_other, goal)
end
