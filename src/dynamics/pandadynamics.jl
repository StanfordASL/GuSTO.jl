export PandaDyn
using PandaRobot

# Useless for now...

mutable struct PandaDyn <: DynamicsModel
  # state: q
  x_dim # joints
  u_dim

  # Parameters that can be updated
  B
end
function PandaDyn()
  x_dim = 9 # 9 joints
  u_dim = 9

  PandaDyn(x_dim,u_dim,[])
end

######
# Cost 
######


##################
# Initializations
##################