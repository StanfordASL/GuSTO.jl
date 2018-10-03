using Quaternions: imag

# Augmenting the Quaternions.jl package
function Base.Vector(q::Quaternion, scalar_first::Bool = true)
	if scalar_first
		return [real(q); imag(q)]
	else
		return [imag(q); real(q)]
	end
end

# Necessary since there is already a competing definition for Quaterion(a::Vector)
# TODO(ambyld): Add a modified version of Quaternion(a::Vector) that does this conversion properly
function vec2quat(qvec::Vector, scalar_first::Bool = true)
	if scalar_first
		return quat(qvec[1], qvec[2:4])
	else
		return quat(qvec[4], qvec[1:3])
	end
end

rotate_vec(v::Vector, q::Quaternion) = imag(q*quat(v)/q)
