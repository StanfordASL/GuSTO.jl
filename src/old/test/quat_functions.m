q_mat = struct('quatmultiply', [], 'quatrotate', [], 'quat2dcm', [],...
                            'dcm2quat', [], 'quatdivide', [])

%% quatmultiply
q_mat.quatmultiply.test1 = quatmultiply([0 1 0 0], [1 0 0 0]);
q_mat.quatmultiply.test2 = quatmultiply([0 1 0 0], [cos(pi/4) sin(pi/4) 0 0]);
q_mat.quatmultiply.test3 = quatmultiply([0 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)], [cos(pi/4) sin(pi/4) 0 0]);

%% quatrotate
q_mat.quatrotate.test1 = quatrotate([1/sqrt(2) -1/sqrt(2) 0 0], [1 1 1]);
q_mat.quatrotate.test2 = quatrotate([-1/sqrt(2) 1/sqrt(2) 0 0], [1 1 1]);
q_mat.quatrotate.test3 = quatrotate([cos(0.25*pi/2) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3)], [1 0 0]);
q_mat.quatrotate.test4 = quatrotate([0.595944  -0.0266435  0.202378 -0.776649], [1 1 1]);

%% quat2dcm
q_mat.quat2dcm.test1 = quat2dcm([1 0 0 0]);
q_mat.quat2dcm.test2 = quat2dcm([0 0 0 1]);
q_mat.quat2dcm.test3 = quat2dcm([-1/sqrt(2) 1/sqrt(2) 0 0]);
q_mat.quat2dcm.test4 = quat2dcm([cos(0.25*pi/2) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3)]);

%% dcm2quat
q_mat.dcm2quat.test1 = dcm2quat(eye(3));
% qstruct.dcm2quat.test2 = dcm2quat(rotx(90)*roty(0)*rotz(90))
% qstruct.dcm2quat.test3 = dcm2quat(rotx(45)*roty(45)*rotz(45))
%qstruct.dcm2quat.test4 = dcm2quat(rotx(45)*roty(0)*rotz(45));
tempq = dcm2quat([0.360941986733097, 0.911127899745766, 0.198914133530110;...
    0.927979080186442, -0.372073016144099, 0.020408267779449;...
    -0.092605123775577, -0.177221953951259, 0.979804403995108]);
q_mat.dcm2quat.test5 = quatnormalize(tempq);

%% quatdivide
q_mat.quatdivide.test1 = quatdivide([1 0 0 0], [1 0 0 0]);
q_mat.quatdivide.test2 = quatdivide([1 0 0 0], [-1 0 0 0]);
q_mat.quatdivide.test3 = quatdivide([cos(0.25*pi/2) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3)], [1 0 0 0]);
q_mat.quatdivide.test4 = quatdivide([0.595944  -0.0266435  0.202378 -0.776649], [cos(0.25*pi/2) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3) sin(0.25*pi/2)/sqrt(3)]);

%% quatinterp
q_mat.quatinterp.test1 = quatinterp([1,0,0,0],[1,0,0,0],0.5,'slerp');
q_mat.quatinterp.test2 = quatinterp([1,0,0,0],[0,-1,0,0],0.5,'slerp');
q_mat.quatinterp.test3 = ...
    quatinterp([0 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)], [cos(pi/4) sin(pi/4) 0 0],0,'slerp');
q_mat.quatinterp.test4 = ...
    quatinterp([0 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)], [cos(pi/4) sin(pi/4) 0 0],1,'slerp');
q_mat.quatinterp.test5 = quatinterp([0 1 0 0], [1 0 0 0], 0.5, 'slerp');

%% quatexp
q_mat.quatexp.test1 = quatexp(0.5*pi*[0,1/sqrt(3),1/sqrt(3),1/sqrt(3)]);
q_mat.quatexp.test2 = quatexp(0.5*pi/3*[0,0,0,-1]);
q_mat.quatexp.test3 = quatexp(0.5*pi/4*[0,0,0,1]);
q_mat.quatexp.test4 = quatexp([0,0,0,0]);

%% quatlog
q_mat.quatlog.test1 = quatlog([cos(pi/3), sin(pi/3)/sqrt(3), sin(pi/3)/sqrt(3), sin(pi/3)/sqrt(3)]);
q_mat.quatlog.test2 = quatlog([cos(pi/2), sin(pi/2)/sqrt(3), sin(pi/2)/sqrt(3), sin(pi/2)/sqrt(3)]);
q_mat.quatlog.test3 = quatlog([cos(pi),0,0,0]);

%% write to file
field_names = fieldnames(q_mat);
for i = 1:numel(field_names)
    field_name = q_mat.(field_names{i});
    tests = fieldnames(field_name);
    for j = 1:numel(tests)
        val = field_name.(tests{j});
        if numel(val) == 4
            % Mat --> Julia quaternion convention
            q_mat.(field_names{i}).(tests{j}) = [val(2:4) val(1)];
        end
    end
end 

save('quat_functions.mat', 'q_mat', 'q_mat');
