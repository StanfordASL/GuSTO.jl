function [RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,qm,robot)
% Computes the kinematics -- positions and orientations -- of the multibody system.
%
% [RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,qm,robot)
%
% :parameters: 
%   * R0 -- Rotation matrix from the base-link CCS to the inertial CCS -- [3x3].
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * qm -- Displacements of the active joints -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%c
% :return: 
%   * RJ -- Joints CCS 3x3 rotation matrices with respect to the inertial CCS  -- as a [3x3xn] matrix.
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * rJ -- Positions of the joints, projected in the inertial CCS -- as a [3xn] matrix.
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * e -- Joint rotation/sliding axes, projected in the inertial CCS -- as a [3xn] matrix.
%   * g -- Vector from the origin of the ith joint CCS to the origin of the ith link CCS, projected in the inertial CCS -- as a [3xn] matrix.
%
% Remember that all the ouput magnitudes are projected in the **inertial frame**.
%
% Examples on how to retrieve the results from a specific link/joint:
%
%   To retrieve the position of the ith link: ``rL(1:3,i)``.
%
%   To retrieve the rotation matrix of the ith joint: ``RJ(1:3,1:3,i)``.   
%
% See also: :func:`src.robot_model.urdf2robot` and :func:`src.robot_model.DH_Serial2robot`.

%{  
    LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

%=== CODE ===%

%--- Number of links and joints ---%
n=robot.n_links_joints;

%--- Homogeneous transformation matrices ---%

%Pre-allocate homogeneous transformations matrices
TJ=zeros(4,4,n,'like',R0);
TL=zeros(4,4,n,'like',R0);

TJ = sym(zeros(4,4,n));
TL = sym(zeros(4,4,n));

%--- Base-link ---%
T0=[R0,r0;zeros(1,3),1];

%--- Forward kinematics recursion ---%

%Obtain the joints and links kinematics
for i=1:n
    
    %Get child joint
    cjoint=robot.joints(i);
    
    %Joint kinematics (homogeneous transformation matrix)
    if cjoint.parent_link==0
        %Parent link is the base-link
        TJ(1:4,1:4,cjoint.id)=T0*cjoint.T;
    else
        %Parent link is not the base-link
        TJ(1:4,1:4,cjoint.id)=TL(1:4,1:4,cjoint.parent_link)*cjoint.T;
    end
    
    %Transformation due to current joint variable
    if cjoint.type==1
        %Revolute
        T_qm=[Euler_DCM(cjoint.axis,qm(cjoint.q_id))',zeros(3,1);zeros(1,3),1];
    elseif cjoint.type==2
        %Prismatic
        T_qm=[eye(3),cjoint.axis*qm(cjoint.q_id);zeros(1,3),1];
    else
        %Fixed
        T_qm=[eye(3),zeros(3,1);zeros(1,3),1];
    end
    
    %Link Kinematics (homogeneous transformation matrix)
    clink=robot.links(cjoint.child_link);
    TL(1:4,1:4,clink.id)=TJ(1:4,1:4,clink.parent_joint)*T_qm*clink.T;
end

%--- Rotation matrices, translation, position and other geometric quantities ---%

%Pre-allocate rotation matrices, translation and positions
RJ=zeros(3,3,n,'like',R0);
RL=zeros(3,3,n,'like',R0);
rJ=zeros(3,n,'like',R0);
rL=zeros(3,n,'like',R0);
%Pre-allocate rotation/sliding axis
e=zeros(3,n,'like',R0);
%Pre-allocate other geometric quantities
g=zeros(3,n,'like',R0);

RJ = sym(zeros(3,3,n));
RL = sym(zeros(3,3,n));
rJ = sym(zeros(3,n));
rL = sym(zeros(3,n));
e = sym(zeros(3,n));
g = sym(zeros(3,n));

%Format rotation matrices, link positions, joint axis and other geometric
%quantities

%Joint associated quantities
for i=1:n
    RJ(1:3,1:3,i)=TJ(1:3,1:3,i);
    rJ(1:3,i)=TJ(1:3,4,i);
    e(1:3,i)=RJ(1:3,1:3,i)*robot.joints(i).axis;
end
%Link associated quantities
for i=1:n
    RL(1:3,1:3,i)=TL(1:3,1:3,i);
    rL(1:3,i)=TL(1:3,4,i);
    g(1:3,i)=rL(1:3,i)-rJ(1:3,robot.links(i).parent_joint);
end


end