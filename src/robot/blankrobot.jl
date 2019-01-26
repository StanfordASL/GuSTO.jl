export BlankRobot

mutable struct BlankRobot <: Robot
  btCollisionObject
end
BlankRobot() = BlankRobot(BulletCollision.convex_hull([zeros(3)]))
