### Unit tests for quat_functions.jl###
### Julia functions are compared against Matlab outputs###

# For each test block:
# First line is Matlab command used to generate value
# Second line is the expected value based on Matlab output (scalar term in Matlab is q(1), here is q[4])
# Third line is Julia command

using MAT
import Base.Test: @test, @testset
include("../quat_functions.jl")

fn = "quat_functions.mat"
file = matopen(fn)
q_mat = read(file, "q_mat")
close(file)




### quat_multiply ###
@testset "quat_multiply" begin
  # Note: Matlab uses reverse order of quaternion multiplication
  
  expected = vec(q_mat["quatmultiply"]["test1"])
  test = quat_multiply([0.;0.;0.;1.], [1.;0.;0.;0.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatmultiply"]["test2"])
  test = quat_multiply([sin(pi/4); 0.; 0.; cos(pi/4)], [1.;0.;0.;0.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatmultiply"]["test3"])
  test = quat_multiply([sin(pi/4); 0; 0; cos(pi/4)], [1/sqrt(3); 1/sqrt(3); 1/sqrt(3); 0])
  @test isapprox(expected,test)
end




### quat_rotate ###
@testset "quat_rotate" begin
  expected = vec(q_mat["quatrotate"]["test1"])
  test = quat_rotate([-1/sqrt(2); 0.; 0.; 1/sqrt(2)], [1.;1.;1.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatrotate"]["test2"])
  test = quat_rotate([1/sqrt(2); 0.; 0.; -1/sqrt(2)], [1.;1;1])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatrotate"]["test3"])
  test = quat_rotate([sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); cos(0.25*pi/2)], [1.; 0.; 0.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatrotate"]["test4"])
  test = quat_rotate([-0.0266435;  0.202378; -0.776649;  0.595944], [1.; 1.; 1.])
  @test isapprox(expected,test)
end




### quat2dcm ###
@testset "quat2dcm" begin
  expected = q_mat["quat2dcm"]["test1"]
  test = quat2dcm([0.; 0.; 0.; 1.])
  @test isapprox(expected,test)
  
  expected = q_mat["quat2dcm"]["test2"]
  test = quat2dcm([0.; 0.; 1.; 0.])
  @test isapprox(expected,test)
  
  expected = q_mat["quat2dcm"]["test3"]
  test = quat2dcm([1/sqrt(2); 0.; 0.; -1/sqrt(2)])
  @test isapprox(expected,test)
  
  expected = q_mat["quat2dcm"]["test4"]
  test = quat2dcm([sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); cos(0.25*pi/2)])
  @test isapprox(expected,test)
end




### dcm2quat ###
@testset "dcm2quat" begin
  # dcm2quat(eye(3))
  expected = [0.; 0.; 0.; 1.]
  expected = vec(q_mat["dcm2quat"]["test1"])
  test = dcm2quat(eye(3))
  @test isapprox(expected,test)
  
  # dcm2quat(rotx(90)*roty(0)*rotz(90))
  # NOTE: Matlab returns negative of expected given here
  expected = -[0.5; -0.5; 0.5; -0.5]
  # expected = vec(q_mat["dcm2quat"]["test2"])
  test = dcm2quat([0. -1. 0.; 0. 0. -1.; 1. 0. 0.])
  @test isapprox(expected,test)
  
  # dcm2quat(rotx(45)*roty(45)*rotz(45))
  expected = [-0.461939766255643; -0.191341716182545; -0.461939766255643; 0.732537816328742]
  # expected = vec(q_mat["dcm2quat"]["test3"])
  test = dcm2quat([0.5 -0.5 0.707106781186547; 0.853553390593274 0.146446609406726  -0.5; 0.146446609406726 0.853553390593274 0.5])
  @test isapprox(expected,test)
  
  # dcm2quat(rotx(45)*roty(0)*rotz(45))
  expected = [-0.353553390593274; 0.146446609406726; -0.353553390593274; 0.853553390593274]
  # expected = vec(q_mat["dcm2quat"]["test4"])
  test = dcm2quat([0.707106781186547 -0.707106781186547 0.; 0.5 0.5 -0.707106781186547; 0.5 0.5 0.707106781186547])
  @test isapprox(expected,test)
  
  # q = dcm2quat([0.360941986733097, 0.911127899745766, 0.198914133530110; 0.927979080186442, -0.372073016144099, 0.020408267779449; -0.092605123775577, -0.177221953951259, 0.979804403995108]);
  # q = q./norm(q);
  # NOTE: Matlab does not return normalized quaternion here
  expected = [0.098815110865354; -0.145759628652843; -0.008425590220338; 0.984336687292052]
  expected = vec(q_mat["dcm2quat"]["test5"])
  test = dcm2quat([0.360941986733097 0.911127899745766 0.198914133530110; 0.927979080186442 -0.372073016144099 0.020408267779449; -0.092605123775577 -0.177221953951259 0.979804403995108])
  @test isapprox(expected,test)
end  




### quat_error ###
# NOTE: in Matlab, d = quatdivide(q,r) s.t. q = quatmultiply(r,d)
@testset "quat_error" begin
  expected = vec(q_mat["quatdivide"]["test1"])
  test = quat_error([0.; 0.; 0.; 1.], [0.; 0.; 0.; 1.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatdivide"]["test2"])
  test = quat_error([0.; 0.; 0.; -1.], [0.; 0.; 0.; 1.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatdivide"]["test3"])
  test = quat_error([0.; 0.; 0.; 1.], [sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); cos(0.25*pi/2)])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatdivide"]["test4"])
  test = quat_error([sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); sin(0.25*pi/2)/sqrt(3); cos(0.25*pi/2)], [-0.0266435;  0.202378; -0.776649; 0.595944])
  @test isapprox(expected,test)
end





### quat_interp ###
@testset "quat_interp" begin
  expected = vec(q_mat["quatinterp"]["test1"])
  test = quat_interp([0; 0; 0; 1.], [0; 0; 0; 1.], 0.5)
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatinterp"]["test2"])
  test = quat_interp([0.; 0; 0; 1], [-1.; 0; 0; 0], 0.5)
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatinterp"]["test3"])
  test = quat_interp([1/sqrt(3); 1/sqrt(3); 1/sqrt(3); 0], [sin(pi/4); 0; 0; cos(pi/4)], 0.) 
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatinterp"]["test4"])
  test = quat_interp([1/sqrt(3); 1/sqrt(3); 1/sqrt(3); 0], [sin(pi/4); 0; 0; cos(pi/4)], 1.) 
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatinterp"]["test5"])
  test = quat_interp([1.; 0; 0; 0], [0.; 0; 0; 1], 0.5)
  @test isapprox(expected,test)
end




### quat_exp ###
@testset "quat_exp" begin
  # NOTE: Matlab appends 0 in element 1 of input vector 
  # NOTE: Matlab uses theta/2*v as input instead of theta*v 
  
  expected = vec(q_mat["quatexp"]["test1"])
  test = quat_exp(pi*[1/sqrt(3);1/sqrt(3);1/sqrt(3)])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatexp"]["test2"])
  test = quat_exp(pi/3*[0;0;-1.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatexp"]["test3"])
  test = quat_exp(pi/4*[0;0;1.])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatexp"]["test4"])
  test = quat_exp([0;0;0.])
  @test isapprox(expected,test)
end





### quat_log ###
@testset "quat_log" begin
  # NOTE: Matlab uses [cos(th);sin(th)*v] instead of th/2 
  expected = vec(q_mat["quatlog"]["test1"])[1:3]
  test = quat_log([sin(pi/6)/sqrt(3); sin(pi/6)/sqrt(3); sin(pi/6)/sqrt(3); cos(pi/6)])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatlog"]["test2"])[1:3]
  test = quat_log([sin(pi/4)/sqrt(3); sin(pi/4)/sqrt(3); sin(pi/4)/sqrt(3); cos(pi/4)])
  @test isapprox(expected,test)
  
  expected = vec(q_mat["quatlog"]["test3"])[1:3]
  test = quat_log([0;0;0;cos(pi/2)])
  @test isapprox(expected,test)
end
