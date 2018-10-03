export Car

mutable struct Car <: Robot 
  btCollisionObject
end
Car() = Car(BulletCollision.convex_hull([zeros(3)]))
