include("scp/scp_gusto.jl")
include("scp/scp_trajopt.jl")
# include("scp/scp_mao.jl")

function add_constraint_category!(SCPCCat, func, dimtype::Symbol, ind_time, ind_other...)
	func_name = Symbol(split(string(func),".")[end])
	cc = ConstraintCategory(func, dimtype, ind_time, ind_other)
	haskey(SCPCCat, func_name) ? push!(SCPCCat[func_name], cc) : SCPCCat[func_name] = [cc]
end

function add_constraint_category!(SCPCCat, func, goal::G, dimtype::Symbol, ind_other...) where G <: Goal{T} where T
	func_name = Symbol(split(string(func),".")[end])
	cc = ConstraintCategory(func, dimtype, goal.ind_coordinates, ind_other, goal)
	haskey(SCPCCat, func_name) ? push!(SCPCCat[func_name], cc) : SCPCCat[func_name] = [cc]
end
