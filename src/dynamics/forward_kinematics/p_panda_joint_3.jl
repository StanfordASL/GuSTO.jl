export p_panda_joint_3
function p_panda_joint_3(x)
	x1  = x[1]
	x2  = x[2]
	x3  = x[3]
	x4  = x[4]
	x5  = x[5]
	x6  = x[6]
	x7  = x[7]
	x8  = x[8]
	x9  = x[9]
	x10 = x[10]
	x11 = x[11]
	x12 = x[12]
	x13 = x[13]
	x14 = x[14]

	A0 = zeros(3)
  A0[0+1,1+0] = x1*x4*(7.9E1/2.5E2)+x2*x3*1.547322079806372E-12;
  A0[1+1,1+0] = x1*x3*(-1.547322079806372E-12)+x2*x4*(7.9E1/2.5E2);
  A0[2+1,1+0] = x3*(7.9E1/2.5E2)+3.33E2/1.0E3;


	return A0
end
