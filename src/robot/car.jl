export Car

mutable struct Car <: Robot 
  btCollisionObject
end
Car() = Car(BulletCollision.sphere(SVector{3}(zeros(Float32,3)), 0.001))
# Necessary - Avoids zero distance when calling Bullet