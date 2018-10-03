export BlankRobot

mutable struct BlankRobot <: Robot
  btCollisionObject
end
BlankRobot() = BlankRobot(BT.convex_hull([zeros(3)]))
