using GeometryTypes
include("../collision_checker.jl")

### TEST 1 ###
TEST_1_passes = true
zones = HyperRectangle[]
corner1,corner2 = [0.1;0.1;0.1], [0.4;0.4;0.4]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
corner1,corner2 = [0.4;0.6;0.0], [0.9;0.9;0.9]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
test_pts = [0 0.25 0.4 0.6 0.6;
            0 0.25 0.4 0.55 0.8;
            0 0.25 0.4 0.4 0.4]
in_free_state = Bool[true,false,false,true,false]

for (k,truth) in enumerate(in_free_state)
  TEST_1_passes &= (truth == is_free(test_pts[:,k],zones))
  !TEST_1_passes && println("Test 1.$(k) fails")
end

TEST_1_passes && println("Test 1 passes")




### TEST 2 ###
TEST_2_passes = true
zones = HyperRectangle[]
corner1,corner2 = [0;0.45;0],[0.9;0.55;0.01]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
corner1,corner2 = [0.25;0.2;0],[0.3;0.8;0.01]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
corner1,corner2 = [0.75;0.2;0],[0.8;0.8;0.01]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
corner1,corner2 = [0.5;0.7;0],[0.55;1;0.01]
push!(zones, HyperRectangle(Vec3f0(corner1), Vec3f0(corner2-corner1)))
corner1,corner2 = [0.5;0;0],[0.55;0.3;0.01]

test_pts = [0.1  0.275  0.275  0.275;
            0.65  0.8  0.7  0.7;
            0  0  0  0.02];
in_free_state = Bool[true,false,false,true]

for (k,truth) in enumerate(in_free_state)
  TEST_2_passes &= (truth == is_free(test_pts[:,k],zones))
  !TEST_2_passes && println("Test 2.$(k) fails")
end

TEST_2_passes && println("Test 2 passes")
