import Base.Test: @test, @testset
import StaticArrays: SVector, MVector
using Cxx
using Revise
using GeometryTypes
include("../bullet_collision.jl")

@testset "gjk_test" begin
  c1 = BT.convex_hull(SVector{3}(-1.,0,0))
  c2 = BT.convex_hull(SVector{3}(4.,0,0))
  @test isapprox(BT.distance(c1,c2)[1],5)

  c1 = BT.convex_hull([SVector{3}(0.,0.,0.), SVector{3}(0.,1.,0.), SVector{3}(1.,0.,0.), SVector{3}(1.,1.,0.)])
  c2 = BT.convex_hull([SVector{3}(4.,0.5,0.)])
  @test isapprox(BT.distance(c1,c2)[1], 3)

  c1 = BT.convex_hull([SVector{3}(1.,2,3)])
  c2 = BT.convex_hull([SVector{3}(3.,4,5)])
  @test isapprox(BT.distance(c1,c2)[1], norm([1.,2,3]-[3.,4,5]))

  c1 = BT.convex_hull([SVector{3}(0,1,0.), SVector{3}(1.,1,0), SVector{3}(0,2.,0),SVector{3}(-1,1.,0)])
  c2 = BT.convex_hull([SVector{3}(0.,0,0)])
  @test isapprox(BT.distance(c1,c2)[1], 1.)
end





@testset "box" begin
  btworld_ = BT.collision_world(SVector{3}(-1000.,-1000.,-1000.), SVector{3}(1000.,1000.,1000.))
  btrobot_ = BT.box(SVector{3}(2.,0.,0.),SVector{3}(3.,1,0))

  btenvironment_ = BT.BulletStaticEnvironment(btrobot_,btworld_)

  BT.add_collision_object!(btenvironment_,BT.convex_hull_box(SVector{3}(0.,0,0),SVector{3}(1.,1,0)))
  @test BT.is_free_state(SVector{3}(2.5,0.5,0.), btenvironment_) == true
  @test BT.is_free_state(SVector{3}(1.51,0.5,0.), btenvironment_) == true # seem to have 1cm margin
  @test BT.is_free_state(SVector{3}(1.52,0.5,0.), btenvironment_) == true
  @test BT.is_free_state(SVector{3}(1.5,0.5,0.), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(1.,0.5,0.), btenvironment_) == false 
end





@testset "sphere" begin
  btworld_ = BT.collision_world(SVector{3}(-1000.,-1000.,-1000.), SVector{3}(1000.,1000.,1000.))
  #btrobot_ = BT.sphere(SVector{3}(0.,0.,0.), 0.3)
  btrobot_ = BT.convex_hull(SVector{3}(0.,0.,0.))
  BT.set_margin(btrobot_,0.3)
  btenvironment_ = BT.BulletStaticEnvironment(btrobot_,btworld_)

  BT.add_collision_object!(btenvironment_,BT.convex_hull_box(SVector{3}(0.,0,0),SVector{3}(1.,1,1)))
  @test BT.is_free_state(SVector{3}(0.5,0.5,0.5), btenvironment_) == false
  @test BT.is_free_state(SVector{3}(1.2,0.5,0.5), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(1.37,0.5,0.5), btenvironment_) == true # seem to have 2cm margin
  @test BT.is_free_state(SVector{3}(1.38,0.5,0.5), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(0.,0.,1.37), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(0.,0.,1.38), btenvironment_) == true 
end





@testset "convex_hull" begin
  btworld_ = BT.collision_world(SVector{3}(-1000.,-1000.,-1000.), SVector{3}(1000.,1000.,1000.))
  btrobot_ = BT.convex_hull_box(SVector{3}(2.,0.,0.),SVector{3}(3.,1,0))
  btenvironment_ = BT.BulletStaticEnvironment(btrobot_,btworld_)

  BT.add_collision_object!(btenvironment_,BT.convex_hull_box(SVector{3}(0.,0,0),SVector{3}(1.,1,0)))
  @test BT.is_free_state(SVector{3}(1.,0.,0.), btenvironment_) == true # seems being on boundary not in collision
  @test BT.is_free_state(SVector{3}(-0.9,0.,0), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(-1.1,0.,0), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(-3.1,0.,0), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(-2.,1.,0), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(-2.,1.01,0), btenvironment_) == false 
  #@test BT.is_free_state(SVector{3}(-2.,1.02,0), btenvironment_) == true 
end





@testset "point" begin
  btworld_ = BT.collision_world(SVector{3}(-1000.,-1000.,-1000.), SVector{3}(1000.,1000.,1000.))
  btrobot_ = BT.convex_hull(SVector{3}(0.,0.,0.))
  btenvironment_ = BT.BulletStaticEnvironment(btrobot_,btworld_)

  BT.add_collision_object!(btenvironment_,BT.convex_hull_box(SVector{3}(0.,0,0),SVector{3}(1.,1,1.)))
  @test BT.is_free_state(SVector{3}(0.,0.,0.), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(0.,1.,0.), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(1.,1.,1.), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(0.5,0.5,0.5), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(0.99,0.01,0.01), btenvironment_) == false 
  @test BT.is_free_state(SVector{3}(1.,0.01,0.01), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(-0.01,0.01,0.01), btenvironment_) == true 
  @test BT.is_free_state(SVector{3}(1.01,0.01,0.01), btenvironment_) == true 
end
