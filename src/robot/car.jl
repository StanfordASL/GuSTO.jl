export Car

mutable struct Car <: Robot 
  btCollisionObject
end
# Necessary - Avoids zero distance when calling Bullet
Car() = Car(BulletCollision.sphere(SVector{3}(zeros(Float32,3)), 0.001))
