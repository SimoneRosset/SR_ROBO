% AUTHOR: SIMONE ROSSETTI
classdef SR_ROBO
    % ROBOTICS TOOLBOX
    properties( Constant = true )
        epsilon = 1e-4;
        epsilon_q = 1e-5;
        min_det = 1e-6;
        max_iter = 50;
    end
    methods(Static)
        %%%%%%%%% ALGEBRA %%%%%%%%%
        
        function help_angle()
            disp('sin(alpha+beta)  =  sin(alpha)*cos(beta)+sin(beta)*cos(alpha);')
            disp('sin(alpha-beta)  =  sin(alpha)*cos(beta)-sin(beta)*cos(alpha);')
            disp('sin(alpha-beta)  =  sin(alpha)*cos(beta)-sin(beta)*cos(alpha);')
            disp('cos(alpha+beta)  =  cos(alpha)*cos(beta)-sin(alpha)*sin(beta);')
            disp('cos(alpha-beta)  =  cos(alpha)*cos(beta)+sin(alpha)*sin(beta);')
            disp('tan(alpha+beta)  =  (tan(alpha)+tan(beta))/(1-tan(alpha)*tan(beta));')
            disp('tan(alpha-beta)  =  (tan(alpha)-tan(beta))/(1+tan(alpha)*tan(beta));')
        end
        
        function m = R_x(gamma)
            % return elementary rotation around X axis
            m =[1       0           0;
                0   cos(gamma)  -sin(gamma);
                0   sin(gamma)  cos(gamma)];
        end
        
        function m = R_y(beta)
            % return elementary rotation around Y axis
            m =[cos(beta)   0    sin(beta);
                0       1       0;
                -sin(beta)  0    cos(beta)];
        end
        
        function m = R_z(alpha)
            % return elementary rotation around Z axis
            m=[cos(alpha)   -sin(alpha)     0;
                sin(alpha)   cos(alpha)      0;
                0           0           1];
        end
        
        function m = R_zyz(rz,ry,rzz)
            % minimal representation - moving axis
            m = SR_ROBO.R_z(rz) * SR_ROBO.R_y(ry) * SR_ROBO.R_z(rzz);
        end
        
        function m = R_zxz(rz,rx,rzz)
            % minimal representation - moving axis
            m = SR_ROBO.R_z(rz) * SR_ROBO.R_x(rx) * SR_ROBO.R_z(rzz);
        end
        
        function m = R_xyz(rx,ry,rz)
            % minimal representation - moving axis
            m = SR_ROBO.R_x(rx) * SR_ROBO.R_y(ry) * SR_ROBO.R_z(rz);
        end
        
        function m = R_xzy(rx,rz,ry)
            % minimal representation - moving axis
            m = SR_ROBO.R_x(rx) * SR_ROBO.R_z(rz) * SR_ROBO.R_y(ry);
        end
        
        function m = R_zyx(rz,ry,rx)
            % minimal representation - moving axis
            m = SR_ROBO.R_z(rz) * SR_ROBO.R_y(ry) * SR_ROBO.R_x(rx);
        end
        
        function m = R_zxy(rz,rx,ry)
            % minimal representation - moving axis
            m = SR_ROBO.R_z(rz) * SR_ROBO.R_x(rx) * SR_ROBO.R_y(ry);
        end
        
        function m = R_yxy(ry,rx,ryy)
            % minimal representation - moving axis
            m = SR_ROBO.R_y(ry) * SR_ROBO.R_x(rx) * SR_ROBO.R_y(ryy);
        end
        
        function m = R_yxz(ry,rx,rz)
            % minimal representation - moving axis
            m = SR_ROBO.R_y(ry) * SR_ROBO.R_x(rx) * SR_ROBO.R_z(rz);
        end
        
        function m = R_rpy(rr,rp,ryy)
            % minimal representation - fixed axis
            [rx,ry,rz] = deal(rr,rp,ryy);
            m = SR_ROBO.R_zyx(rz,ry,rx);
        end
        
        function [alpha,beta,gamma] = R_zyz_inv(R,positive)
            assert(size(R,1)==3 && size(R,2)==3,'R must be 3 x 3')
            assert(R(1,3) ~= 0 && R(2,3) ~= 0, 'singular tern, beta = 0 or beta = pi')
            assert(positive==0 || positive==1, 'you must choose which configuration {1,0} given by the sign of the square root of r13 and r23 components')
            if positive
                alpha = atan2(R(2,3),R(1,3));
                beta = atan2(sqrt(R(2,3)^2+R(1,3)^2),R(3,3));
                gamma = atan2(R(3,2),-R(3,1));
            else
                alpha = atan2(-R(2,3),-R(1,3));
                beta = atan2(-sqrt(R(2,3)^2+R(1,3)^2),R(3,3));
                gamma = atan2(-R(3,2),R(3,1));
            end
        end
        
        function R = ax_ang(theta, axis)
            assert(size(theta,1)==1 && size(theta,2)==1, 'theta must be scalar')
            assert(size(axis,1)==3 && size(axis,2)==1, 'axis must be column vector 3 x 1')
            R = axis * axis.' + (eye(3) - axis * axis.') * cos(theta) + [0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0 ]*sin(theta);
        end
        
        function [angle,axis]= ax_ang_inv(R)
            assert(size(R,1)==3 && size(R,2)==3, 'rotation matrix must be 3 x 3')
            angle = atan2(sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2),(R(1,1)+R(2,2)+R(3,3)-1));
            if sin(angle) ~= 0
                axis = 1/(2*sin(angle)).*[R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]';
            else
                if abs(angle) == pi
                    axis = [sqrt((R(1,1)+1)/2),sqrt((R(2,2)+1)/2),sqrt((R(3,3)+1)/2)]';
                else
                    axis = [0 0 0]'; % undefined
                end
            end
        end
        %%%%%%%%% DK %%%%%%%%%
        function showDHTable(matrix)
            % display the DH table
            cell_=cell(size(matrix,1),size(matrix,2)+1);
            for i=1:size(matrix,1)
                cell_{i,1}=i;
                for j=2:size(matrix,2)+1
                    cell_{i,j}=char(matrix(i,j-1));
                end
            end
            disp(cell2table(cell_,...
                'VariableNames',{'joint_i','alpha_i' 'a_i' 'd_i' 'theta_i'}))
        end
        function draw(A_list,show, varargin)
            % MANIPULATOR DK DRAWING
            % Input: list of T matrix for each joint
            % Output: plot in 3D and compute DK of a robot
            % DK = [N,S,A|P]
            T0 = eye(4,4);
            if length(varargin)==2 && isequal(varargin{1},'T_reference') && isequal(size(T0),[4,4])
                T0 = varargin{2};
                disp('BASE FRAME: ');
                disp(T0)
            end
            DK = T0;
            DK(1:3,1:3) = DK(1:3,1:3)*1/4;
            quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,1),DK(2,1),DK(3,1),'r','LineWidth',2)
            hold on
            quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,2),DK(2,2),DK(3,2),'g','LineWidth',2)
            quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,3),DK(2,3),DK(3,3),'b','LineWidth',2)
            if show
                text(DK(1,4),DK(2,4),DK(3,4), 'BASE')
                text(DK(1,1),DK(2,1),DK(3,1), 'x0','Color',[1, 0 ,0])
                text(DK(1,2),DK(2,2),DK(3,2), 'y0','Color',[0, 1 ,0])
                text(DK(1,3),DK(2,3),DK(3,3), 'z0','Color',[0, 0 ,1])
            end
            DK(1:3,1:3) = DK(1:3,1:3)*4;
            min_=0;
            max_=0;
            for i=1:length(A_list)-1
                DKprev=DK;
                DK = DK * A_list{i};
                DK(1:3,1:3) = DK(1:3,1:3)*1/4;
                quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,1),DK(2,1),DK(3,1),'r','LineWidth',2)
                quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,2),DK(2,2),DK(3,2),'g','LineWidth',2)
                quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,3),DK(2,3),DK(3,3),'b','LineWidth',2)
                if show
                    text(DK(1,4),DK(2,4),DK(3,4), ['link', num2str(i)])
                    text(DK(1,4)+DK(1,1),DK(2,4)+DK(2,1),DK(3,4)+DK(3,1), ['x',num2str(i)],'Color',[1, 0 ,0])
                    text(DK(1,4)+DK(1,2),DK(2,4)+DK(2,2),DK(3,4)+DK(3,2), ['y',num2str(i)],'Color',[0, 1 ,0])
                    text(DK(1,4)+DK(1,3),DK(2,4)+DK(2,3),DK(3,4)+DK(3,3), ['z',num2str(i)],'Color',[0, 0 ,1])
                end
                plot3(DKprev(1,4),DKprev(2,4),DKprev(3,4), 'ko')
                plot3([DKprev(1,4),DK(1,4)],[DKprev(2,4),DK(2,4)],[DKprev(3,4),DK(3,4)], 'k','LineWidth',1)
                DK(1:3,1:3) = DK(1:3,1:3)*4;
                if min(min(DK))<min_
                    min_=min(min(DK));
                end
                if max(max(DK))>max_
                    max_=max(max(DK));
                end
                if show
                    disp(['A', num2str(i),' =']);
                    disp(A_list{i});
                end
            end
            DKprev=DK;
            DK = DK * A_list{length(A_list)};
            DK(1:3,1:3) = DK(1:3,1:3)*1/4;
            h1 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,1),DK(2,1),DK(3,1),'r','LineWidth',2);
            h2 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,2),DK(2,2),DK(3,2),'g','LineWidth',2);
            h3 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,3),DK(2,3),DK(3,3),'b','LineWidth',2);
            if show
                text(DK(1,4),DK(2,4),DK(3,4), 'EE')
                text(DK(1,4)+DK(1,1),DK(2,4)+DK(2,1),DK(3,4)+DK(3,1), ['x',num2str(length(A_list))],'Color',[1, 0 ,0])
                text(DK(1,4)+DK(1,2),DK(2,4)+DK(2,2),DK(3,4)+DK(3,2), ['y',num2str(length(A_list))],'Color',[0, 1 ,0])
                text(DK(1,4)+DK(1,3),DK(2,4)+DK(2,3),DK(3,4)+DK(3,3), ['z',num2str(length(A_list))],'Color',[0, 0 ,1])
            end
            plot3(DKprev(1,4),DKprev(2,4),DKprev(3,4), 'ko')
            plot3([DKprev(1,4),DK(1,4)],[DKprev(2,4),DK(2,4)],[DKprev(3,4),DK(3,4)], 'k','LineWidth',1)
            DK(1:3,1:3) = DK(1:3,1:3)*4;
            if min(min(DK))<min_
                min_=min(min(DK));
            end
            if max(max(DK))>max_
                max_=max(max(DK));
            end
            if show
                disp(['A',num2str(length(A_list)),' =']);
                disp(A_list{length(A_list)});
                disp(['T0',num2str(length(A_list)),' =']);
                if ~ isnumeric(DK)
                    DK = simplify(DK);
                end
                disp(DK);
            end
            legend([h1 h2 h3],{'x','y','z'});
            xlim([(min_-1) (max_+1)]);
            ylim([min_-1 max_+1]);
            zlim([min_-1 max_+1]);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            view([110 15]);
            pbaspect([1 1 1]);
            hold off
        end
        
        function workspace(T_sym,old,new)
            % MANIPULATOR WORKSPACE DRAWING
            % Input: symbolical DK matrix T_sym, variable symbols, values of symbols
            % Output: plot in 3D approximation of workspace
            xM = double(subs(T_sym(1,4),old,new));
            yM = double(subs(T_sym(2,4),old,new));
            zM = double(subs(T_sym(3,4),old,new));
            if T_sym(3,4)~=0
                [C,v] = convhull(xM(:),yM(:),zM(:));
                trisurf(C,xM(:),yM(:),zM(:), ...
                    'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.05, "LineStyle","none");
            else
                C = delaunay(xM(:),yM(:));
                trisurf(C,xM(:),yM(:),zM(:), ...
                    'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1, "LineStyle","none");
            end
            hold on
            plot3(0,0,0,'-ro', 'MarkerSize',12, "LineWidth",5);
            %plot3(xM(:),yM(:),zM(:),'bo');
            hold off
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            view([110 15]);
            pbaspect([1 1 1]);
        end
        
        function m = dh2mat(alpha,a,d,theta)
            % convert params to matrix according to DH convention
            m = [cos(theta) -cos(alpha)*sin(theta) sin(alpha)*sin(theta) a*cos(theta);
                sin(theta) cos(alpha)*cos(theta) -sin(alpha)*cos(theta) a*sin(theta);
                0            sin(alpha)            cos(alpha)            d;
                0                0                     0                 1];
        end
        
        function [DK, A_list, type] = dh2dk(DH_table, show)
            % DH PARAMETERS TABLE TO DK LIST CONVERSION: RETURN (AND SHOW) A LIST AND FULL DK
            % Input: DH table
            % Output: full DK and A list, if show=true display the A list and full DK
            DK = eye(4);
            A_list = cell(1,size(DH_table,1));
            type = {};
            for i = 1:size(DH_table,1)
                alpha = DH_table(i,1);
                a = DH_table(i,2);
                d = DH_table(i,3);
                theta = DH_table(i,4);
                A_list{i} = SR_ROBO.dh2mat(alpha,a,d,theta);
                if theta==0
                    type{i}='p';
                else
                    type{i}='r';
                end
                if show
                    disp(['A', num2str(i),' =']);
                    disp(A_list{i});
                end
                DK = DK * A_list{i};
            end
            if show
                disp(['T0',num2str(length(A_list)),' =']);
                if ~ isnumeric(DK)
                    DK = simplify(DK);
                end
                disp(DK);
            end
        end
        
        %%%%%%%%% IK %%%%%%%%%
        function q_opt = grad(F_sym,q0,rd,alpha)
            % GRADIENT MAX DESCENT METHOD FOR INVERSE KINEMATIC
            % Input: simbolical DK matrix F_sym, initial angles q0, desired
            % target rd, scalar step size alpha > 0
            % Output: q_opt optimal joint values
            assert(isa(alpha,'double'),'alpha is not type double.')
            assert(alpha>0,'alpha should be greater than 0.')
            iter = 0;
            J_q = jacobian(F_sym,symvar(F_sym));
            e = (rd - F_sym);
            J_k = double(subs(J_q,sym2cell(symvar(J_q)),num2cell(q0')));
            e_k = double(subs(e,sym2cell(symvar(e)),num2cell(q0')));
            J_e = J_k' * e_k;
            if sum(abs(J_e),'all') < SR_ROBO.min_det && norm(e_k) > SR_ROBO.epsilon
                disp('e is in the null space of J');
                fprintf('iter : %d \n', iter);
                q_opt = q0;
                return
            end
            q_k = q0;
            q_kk = q0 + alpha * J_e; % get stuck if e_k is in null space of J_k
            while (norm(q_kk - q_k) > SR_ROBO.epsilon_q &&  norm(e_k) > SR_ROBO.epsilon)
                if iter >= SR_ROBO.max_iter
                    disp('max iter exceeded');
                    break
                end
                q_k = q_kk;
                J_k = double(subs(J_q,sym2cell(symvar(J_q)),num2cell(q_k')));
                e_k = double(subs(e,sym2cell(symvar(e)),num2cell(q_k')));
                J_e = J_k' * e_k;
                if sum(abs(J_e))< SR_ROBO.min_det && norm(e_k) > SR_ROBO.epsilon
                    disp('e is in the null space of J');
                    break
                end
                q_kk = q_k + alpha * J_e;
                iter = iter + 1;
            end
            fprintf('iter : %d \n', iter);
            fprintf('error : %f \n', norm(e_k));
            disp('q_opt : ');
            disp(q_kk);
            q_opt = q_kk;
        end
        function q_opt = newton(F_sym,q0,rd)
            % NEWTON METHOD FOR INVERSE KINEMATIC
            % Input: simbolical DK matrix F_sym, initial angles q0, desired
            % target rd
            % Output: q_opt optimal joint values
            iter = 0;
            J = jacobian(F_sym,symvar(F_sym));
            [rows, cols] = size(J);
            if rank(double(subs(J,sym2cell(symvar(J)),num2cell(q0')))) ~= min([rows cols])
                disp('J is singular');
                fprintf('iter : %d \n', iter);
                q_opt = q0;
                return
            end
            e = (rd - F_sym);
            J_inv = inv(J);
            J_k = double(subs(J_inv,sym2cell(symvar(J)),num2cell(q0')));
            e_k = double(subs(e,sym2cell(symvar(e)),num2cell(q0')));
            J_e = J_k * e_k;
            if sum(abs(J_e),'all')< SR_ROBO.min_det && norm(e_k) > SR_ROBO.epsilon
                disp('e is in the null space of J');
                fprintf('iter : %d \n', iter);
                q_opt = q0;
                return
            end
            q_k = q0;
            q_kk = q0 + J_e;
            while norm(q_kk - q_k) > SR_ROBO.epsilon_q && norm(e_k) > SR_ROBO.epsilon
                if iter >= SR_ROBO.max_iter
                    disp('max iter exceeded');
                    break
                end
                q_k = q_kk;
                if rank(double(subs(J,sym2cell(symvar(J)),num2cell(q_k')))) ~= min([rows cols])
                    disp('J is singular');
                    break
                end
                J_k = double(subs(J_inv,sym2cell(symvar(J)),num2cell(q_k')));
                e_k = double(subs(e,sym2cell(symvar(e)),num2cell(q_k')));
                J_e = J_k * e_k;
                if sum(abs(J_e),'all')< SR_ROBO.min_det && norm(e_k) > SR_ROBO.epsilon
                    disp('e is in the null space of J');
                    break
                end
                q_kk = q_k + J_e;
                iter = iter + 1;
            end
            fprintf('iter : %d \n', iter);
            fprintf('error : %f \n', norm(e_k));
            disp('q_opt : ');
            disp(q_kk);
            q_opt = q_kk;
        end
        function twojointspace(F_sym,rd)
            % PLOT SOLUTION SPACE OF 2 LINK ARM
            % Input: simbolical DK matrix F_sym, desired target rd
            % Output: solution space plot
            e = (rd - F_sym);
            [o1,o2]=meshgrid(-pi:0.5:pi);
            e_ = e'*e;
            er = double(subs(e_,sym2cell(symvar(e_)),{o1, o2}));
            contour3(er,200);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            view([110 15]);
            pbaspect([1 1 1]);
        end
        
        %%%%%%%%% DIFF K %%%%%%%%%
        function s = skew(dx,dy,dz)
            % return skew matrix
            s = [0      -dz     dy;
                dz      0       -dx;
                -dy     dx      0];
        end
        
        function bool = isSkew(S)
            % return skew matrix and display omega
            bool = 1-sum(S+S','all');
        end
        
        function w = getWfromSkew(S)
            % return skew matrix and display omega
            if SR_ROBO.isSkew(S)
                w = [-S(2,3) S(1,3) -S(1,2)].';
            else
                w = false;
            end
        end
        
        function [T,angles] = getMappingT(angle1,angle2,angle3,sequence)
            % takes moving axis notation
            % w = T * [alpha_dot,beta_dot,gamma_dot] return T
            assert(isa(angle1,'double') || isa(angle1,'sym'),'r1 is not double nor symbolic');
            assert(isa(angle2,'double') || isa(angle2,'sym'),'r2 is not double nor symbolic');
            assert(isa(angle3,'double') || isa(angle3,'sym'),'r3 is not double nor symbolic');
            assert(strlength(sequence)==3,'rotations sequence in 3D can be only 3');
            sequence=(split(lower(sequence),''));
            sequence=sequence(~cellfun('isempty',sequence));
            sequence=char(sequence);
            if ~isequal(sequence,['r','p','y']')
                assert(sequence(1)=='x' || sequence(1)=='y' || sequence(1)=='z','sequence can contain only euler characters x, y, z, or only sequence rpy');
                assert(sequence(2)=='x' || sequence(2)=='y' || sequence(2)=='z','sequence can contain only euler characters x, y, z, or only sequence rpy');
                assert(sequence(3)=='x' || sequence(3)=='y' || sequence(3)=='z','sequence can contain only euler characters x, y, z, or only sequence rpy');
            end
            x=[1 0 0]';
            y=[0 1 0]';
            z=[0 0 1]';
            R=eye(3,3);
            T=[0 0 0]';
            [r1,r2,~]= deal(angle1,angle2,angle3);
            angles = diff([angle1, angle2, angle3]');
            sequence_temp = sequence;
            if isequal(sequence,['r','p','y']')
                sequence = ['z','y','x']';
                [~,r2,r1]= deal(angle1,angle2,angle3);
            end
            for i=1:3
                if sequence(i)=='x'
                    if i==1
                        T=x;
                        R=SR_ROBO.R_x(r1);
                    elseif i==2
                        T=[T,R*x];
                        R=R*SR_ROBO.R_x(r2);
                    elseif i==3
                        T=[T,R*x];
                    end
                elseif sequence(i)=='y'
                    if i==1
                        T=y;
                        R=SR_ROBO.R_y(r1);
                    elseif i==2
                        T=[T,R*y];
                        R=R*SR_ROBO.R_y(r2);
                    elseif i==3
                        T=[T,R*y];
                    end
                elseif sequence(i)=='z'
                    if i==1
                        T=z;
                        R=SR_ROBO.R_z(r1);
                    elseif i==2
                        T=[T,R*z];
                        R=R*SR_ROBO.R_z(r2);
                    elseif i==3
                        T=[T,R*z];
                    else
                    end
                end
            end
            if isequal(sequence_temp,['r','p','y']')
                temp = T(:,1);
                T(:,1) = T(:,3);
                T(:,3) = temp;
            end
        end
        
        function j = geometric_jacobian(DH_table)
            % compute the geometric jacobian relative to a DH convention
            % manipulator
            [T_sym, A_list_sym, type] = SR_ROBO.dh2dk(DH_table,false);
            j = sym(zeros(6,size(DH_table,1)));
            A = eye(4,4);
            z0 = [0 0 1]';
            p0E = T_sym(1:3,4);
            if strcmp(type{1},'r')
                j(1:3,1)=cross(z0,p0E);
                j(4:6,1)=z0;
            else
                j(1:3,1)=z0;
            end
            for i = 2:size(DH_table,1)
                A = A * A_list_sym{i-1};
                z = A(1:3,1:3) * z0;
                p = p0E - A(1:3,4);
                if strcmp(type{i},'r')
                    j(1:3,i)=cross(z,p);
                    j(4:6,i)=z;
                else
                    j(1:3,i)=z;
                end
            end
            j=simplify(j);
        end
        
        function J_ = jacobian_derivative(J, q)
            % Thanks to Salvatore Cognetta
            % J = Jacobian matrix
            % q = array of sym
            
            syms q_dot_1 q_dot_2 q_dot_3 q_dot_4 q_dot_5 q_dot_6 q_dot_7 real %derivative of qi
            q_ = [q_dot_1 q_dot_2 q_dot_3 q_dot_4 q_dot_5 q_dot_6 q_dot_7]'; %up to 7 uknown joint velocity
            
            nums = size(J);
            J_ = sym(eye(nums));
            
            cell = 0;
            for i = 1:nums(1)
                for k = 1:nums(2)
                    for j = 1:nums(2)
                        cell = cell + diff(J(i,j),q(k))*q_(j);
                    end
                    J_(i,k) = cell;
                    cell = 0;
                end
            end
            J_ = simplify(J_)
        end
        
        %%%%%%%%% DIFF IK %%%%%%%%%
        function draw_clean(A_list_sym,old,new,radius,keep)
            % MANIPULATOR DIFF IK DRAWING
            % Input: list of T matrix for each joint
            % Output: plot in 3D and compute motion plot of a robot
            % DK = [N,S,A|P]
            if keep
                clf
            end
            hold on
            DK = eye(4,4);
            for i=1:length(A_list_sym)-1
                DKprev=DK;
                DK = DK * double(subs(A_list_sym{i},old,new));
                DK(1:3,1:3)=DK(1:3,1:3)*1/4;
                plot3([DKprev(1,4),DK(1,4)],[DKprev(2,4),DK(2,4)],[DKprev(3,4),DK(3,4)], 'k','LineWidth',1)
                if keep
                    quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,1),DK(2,1),DK(3,1),'r','LineWidth',2);
                    quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,2),DK(2,2),DK(3,2),'g','LineWidth',2);
                    quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,3),DK(2,3),DK(3,3),'b','LineWidth',2);
                end
                DK(1:3,1:3) = DK(1:3,1:3)*4;
                
            end
            DKprev=DK;
            DK = DK * double(subs(A_list_sym{length(A_list_sym)},old,new));
            DK(1:3,1:3)=DK(1:3,1:3)*1/4;
            plot3([DKprev(1,4),DK(1,4)],[DKprev(2,4),DK(2,4)],[DKprev(3,4),DK(3,4)], 'k','LineWidth',1)
            h1 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,1),DK(2,1),DK(3,1),'r','LineWidth',2);
            h2 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,2),DK(2,2),DK(3,2),'g','LineWidth',2);
            h3 = quiver3(DK(1,4),DK(2,4),DK(3,4),DK(1,3),DK(2,3),DK(3,3),'b','LineWidth',2);
            xlim([-radius +radius]);
            ylim([-radius +radius]);
            zlim([-radius +radius]);
            legend([h1 h2 h3],{'x','y','z'});
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            grid on
            view([0 90]);
            pbaspect([1 1 1]);
            hold off
        end
        
        function j = dls_jacobian(jacobian, lamb)
            % return damped least square solution matrix
            m = size(jacobian,1);
            j = jacobian'/(lamb*eye(m,m)+jacobian*jacobian');
        end
        
        function j = pseudo_inverse_jacobian(jacobian)
            % return the pseudo inverse of a matrix
            j = jacobian'/(jacobian*jacobian');
        end
        
        function draw_manipulability_ellipsoid(A, r)
            % A should be JJ' for velocity and inv(JJ') for force
            assert(isequal(size(r),[2 1]) || isequal(size(r),[3 1]),'r should be 2D or 3D column vector');
            if size(r,1)<3
                r=[r;0];
            end
            [U,D,V] = svd(A);
            eigvals = diag(D);
            sigma=sqrt(double(eigvals));
            plot3(r(1),r(2),r(3),'r+','LineWidth',1)
            dirx = sigma(1)/2*V(:,1);
            if length(sigma)==2
                sigma=[sigma;0];
                dirx = [dirx; 0];
                dirx(1) = -dirx(1);
            end
            
            [x, y, z] = ellipsoid(0,0,0,sigma(1)/2,sigma(2)/2,sigma(3)/2);
            az = double(atan2d(dirx(2),dirx(1)));
            ay = double(atan2d(dirx(3),dirx(1)));
            ax = double(atan2d(dirx(3),dirx(2)));
            S = surf(x, y, z);
            if abs(az)>1e-5
                rotate(S,[0 0 1],az);
            end
            if abs(ay)>1e-5
                rotate(S,[0 1 0],ay);
            end
            if abs(ax)>1e-5
                rotate(S,[1 0 0],ax);
            end
            S.XData = S.XData+r(1);
            S.YData = S.YData+r(2);
            S.ZData = S.ZData+r(3);
            set(S, 'EdgeColor','none', 'FaceAlpha',0.2);
            
            for i=1:size(V,1)
                dir = sigma(i)/2*V(:,i);
                if size(dir,1)<3
                    dir=[dir;0];
                end
                quiver3(r(1),r(2),r(3),dir(1),dir(2),dir(3),'r','LineWidth',1)
            end
            
        end
        
        %%%%%%%%% T P %%%%%%%%%
        % joint
        function q = ptp_cubic(q0,q1,q_dot_0,q_dot_1,lamb)
            % return PTP cubic polynomial value for a given task and a given
            % timing law, 4 constraints, 3 degree
            % in spline composition (doubly normalized)
            assert(lamb>=0 && lamb<=1, 'timing law parameter lambda should be in range [0,1]')
            Dq = q1-q0;
            c = q_dot_0/Dq;
            d = 0;
            if q_dot_0==q_dot_1 && q_dot_0==0
                a = -2;
                b = 3;
            else
                a = q_dot_0/Dq + q_dot_1/Dq - 2;
                b = 1 - a - c;
            end
            q = q0 + Dq * (a*lamb^3+b*lamb^2+c*lamb+d);
        end
        
        function q = ptp_quintic(q0,q1,q_dot_0,q_dot_1,q_ddot_0,q_ddot_1,lamb)
            % return PTP quintic polynomial value for a given task and a given
            % timing law, 6 constraints, 5 degree
            % in spline composition (doubly normalized)
            assert(lamb>=0 && lamb<=1, 'timing law parameter lambda should be in range [0,1]')
            Dq = q1-q0;
            if q_dot_0==q_dot_1 && q_ddot_0==q_ddot_1 && q_dot_0==q_ddot_0 && q_dot_0==0
                q = q0 + Dq * (6*lamb^5-15*lamb^4+10*lamb^3);
            else
                q = (1-lamb)^3*(q0+(3*q0+q_dot_0)*lamb + ...
                    (q_ddot_0+6*q_dot_0+12*q0)*lamb^2/2) + ...
                    lamb^3*(q1+(3*q1-q_dot_1)*(1-lamb) + ...
                    (q_ddot_1-6*q_dot_1+12*q1)*(1-lamb)^2/2);
            end
        end
        
        function f = rest_to_rest_higher_order(order)
            % compute ptp poly coefficients for higher derivative (initial 
            % and final set to 0) than position for whatever
            % derivative order
            assert(mod(order,1) == 0, 'order must be integer')
            assert(mod(order,2)~=0, 'order must be odd')
            a = flip(sym('a',[order+1 1]));
            l = sym('l','real');
            f = 0;
            for i=1:length(a)
                f=f+a(length(a)+1-i)*l^(i-1);
            end
            A = zeros(order+1,order+1);
            b = zeros(order+1, 1);
            f_0 = subs(f,{l},0);
            row_const_0 = double(gradient(f_0, a)');
            A(1,:) = row_const_0;
            b(1)=0;
            f_1 = subs(f,{l},1);
            row_const_1 = double(gradient(f_1, a)');
            A(2,:)= row_const_1;
            b(2)=1;
            for i=1:fix(order/2)
                der = diff(f,l,i);
                der_0 = subs(der,{l},0);
                row_const_0 = double(gradient(der_0, a)');
                A(i*2+1,:)= row_const_0;
                b(i*2+1)=0;
                der_1 = subs(der,{l},1);
                row_const_1 = double(gradient(der_1, a)');
                A(i*2+2,:)= row_const_1;
                b(i*2+2)=0;
            end
            sol = linsolve(A,b);
            f = subs(f,a,sol);
        end
        
        function A = spline_matrix_A(h)
            % spline coefficient support method
            n = length(h)+1;
            main = zeros(1,n-2);
            for i=1:n-2
                main(i)= 2*(h(i)+h(i+1));
            end
            up = h(1:n-3);
            down = h(3:n-1);
            
            A = diag(main) + [zeros(n-3,1) diag(up); zeros(1,n-2)] ...
                + [zeros(1,n-2); diag(down) zeros(n-3,1)];
        end
        
        function b = spline_constant_b(h,q,v1,vn)
            % spline coefficient support method
            n = length(q);
            b = zeros(n-2,1);
            for i=1:n-2
                b(i)=(3/(h(i)*h(i+1)))* ...
                    (h(i)^2*(q(i+2)-q(i+1))+ ...
                    h(i+1)^2*(q(i+1)-q(i)));
            end
            b(1)=b(1)-h(2)*v1;
            b(n-2)=b(n-2)-h(n-2)*vn;
        end
        
        function a = cubic_spline_coefficients(knots)
            % generate spline coefficients
            assert(length(knots)>2,'knots should be at least 3');
            n = length(knots);
            a = zeros(1,(n-1)*4);
            h = zeros(1,n-1);
            q = zeros(1,n);
            v = zeros(1,n);
            for i=1:n
                if i<n
                    h(i)=knots{i+1}{1}-knots{i}{1};
                end
                q(i)=knots{i}{2};
            end
            v1 = knots{1}{3};
            vn = knots{n}{3};
            A = SR_ROBO.spline_matrix_A(h);
            b = SR_ROBO.spline_constant_b(h,q,v1,vn);
            v(1)=v1;
            v(n)=vn;
            v(2:n-1) = linsolve(A,b);
            for k=1:n-1
                j = 4*k-3;
                a(j) = q(k);
                a(j+1) = v(k);
                C = [h(k)^2 h(k)^3; 2*h(k) 3*h(k)^2];
                d = [q(k+1)-q(k)-v(k)*h(k); v(k+1)-v(k)];
                a(j+2:j+3) = linsolve(C,d);
            end
        end
        
        function q = cubic_spline(coeff,knots,lamb)
            % translate coefficients and time in spline trajectory
            assert(lamb>=0 && lamb<=1, 'timing law parameter lambda should be in range [0,1]');
            assert(length(coeff) == 4*(length(knots)-1), 'knots and condition should multiple by a factor of 4');
            n = length(knots);
            for i=1:n-1
                if lamb>=knots{i}{1} && lamb<=knots{i+1}{1}
                    j = 4*i-3;
                    [a0, a1, a2, a3]=deal(coeff(j),coeff(j+1),coeff(j+2),coeff(j+3));
                    lamb = (lamb - knots{i}{1});
                    break;
                end
            end
            q = a0 + a1 * lamb + a2 * lamb^2 + a3 * lamb^3;
        end
        
        % cartesian
        function p = ptp_linear(p0,p1,s)
            Dp = p1-p0;
            p = p0 + s*Dp;
        end
        
        function [s,T,Ts] = ptp_bcb(L,a_max,v_max,t)
            % perform BANG COAST BANG between two points
            % L: distance between points
            % a_max: max acceleration
            % v_max: max velocity
            % t: current time
            % advice: use this to get to know resulting time interval before
            % [~,T,Ts] = ptp_bcb(L,a_max,v_max,0)
            
            assert(L>0,'L must be positive')
            assert(a_max>0, 'a_max must be positive')
            assert(v_max>0, 'v_max must be positive')
            if L > v_max^2/a_max % there is coast interval -> trapezoid velocity
                Ts = v_max/a_max;
                T = (L*a_max+v_max^2)/(a_max*v_max);
                assert(t<=T,'t must be lower than duration T')
                if t>=0 && t<=Ts
                    s = a_max*t^2*0.5;
                elseif t>=Ts && t<=T-Ts
                    s = v_max*t - v_max^2/a_max * 0.5;
                elseif t>=T-Ts && t<=T
                    s = - a_max*(t-T)^2 * 0.5 + v_max*T - v_max^2/a_max;
                end
            else % no coast -> spike velocity
                T = sqrt(2*L/a_max);
                Ts = T;
                if t>=0 && t<=T/2
                    s = a_max*t^2*0.5;
                else
                    s =  a_max*(T/2)^2 - a_max*(t-T)^2 * 0.5;
                end
            end
        end
        
        function [p,DT] = over_fly(A,B,C,v1,v2,t,varargin)
            % perform different kinds of overfly depending of given
            % conditions, i.e. distance from curvature angle, max
            % acceleration, interval, ..
            assert(length(A)==2 || length(A)==3,'vectors must be 2D or 3D')
            assert(length(A)==length(B)&&length(A)==length(C),'{A,B,C} vectors should have same length')
            assert(isa(v1,'double') && isa(v2,'double'), '{v1,v2} must be scalar')
            assert(isa(t,'double'), '{t} must be scalar')
            assert(length(varargin)==2, 'one condition of {d1,d2,DT,a_max} must be specified')
            
            if strcmp('d1', varargin{1})
                d1 = varargin{2};
                assert(isa(d1,'double'), '{d1} must be scalar')
            elseif strcmp('d2', varargin{1})
                d2 = varargin{2};
                assert(isa(d2,'double'), '{d2} must be scalar')
            elseif strcmp('DT', varargin{1})
                DT = varargin{2};
                assert(isa(DT,'double'), '{DT} must be scalar')
            elseif strcmp('a_max', varargin{1})
                a_max = varargin{2};
                assert(isa(a_max,'double'), '{a_max} must be scalar')
            end
            
            assert(exist('d1','var') || exist('d2','var') || ...
                exist('DT','var') || exist('a_max','var'),...
                'one of {d1,d2,DT,a_max} condition must be specified')
            
            K_ab = (B-A)/norm(B-A);
            K_bc = (C-B)/norm(C-B);
            
            if exist('d1','var')
                DT=2*d1/v1;
                d2=d1*v2/v1;
            elseif exist('d2','var')
                DT=2*d2/v2;
                d1=d2*v1/v2;
            elseif exist('DT','var')
                d1=v1*DT*0.5;
                d2=v2*DT*0.5;
            elseif exist('a_max','var')
                DT = norm(K_bc.*v2-K_ab.*v1)/a_max;
                d1=v1*DT*0.5;
                d2=v2*DT*0.5;
            end
            
            AA = B - K_ab * d1;
            CC = B + K_bc * d2;
            
            p = AA + v1 * K_ab * t + ( v2 * K_bc - v1 * K_ab ) * t^2 / ( 2 * DT );
            
        end
        
        function [p,c,r] = circular_3_points(A,B,C,s)
            % draw a circular path in space given three points
            assert(length(A)==2 || length(A)==3,'vectors must be 2D or 3D')
            assert(length(A)==length(B)&&length(A)==length(C),'{A,B,C} vectors should have same length')
            d = length(A);
            if d==2
                A=[A;0];
                B=[B;0];
                C=[C;0];
            end
            u1= B-A;
            w1=cross((C-A),u1);
            u = u1/norm(u1);
            w = w1/norm(w1);
            v = cross(w,u);
            % u, v orthogonal vectors spanning the plane
            [bx,~]=deal(dot((B-A),u),0);
            [cx,cy]=deal(dot((C-A),u),dot((C-A),v));
            h = ((cx-bx/2)^2+cy^2-(bx/2)^2)/(2*cy);
            c= A+(bx/2)*u+h*v;
            r = norm(B-c);
            if d==2
                c=c(1:2);
            end
            theta = s * 2 * pi;
            p = c + r * cos(theta) * u  + r * sin(theta) * v;
        end
        
        function [t,n,b,k,p_s]= frenet(p,s)
            % p timing law (trajectory)
            % s parameter 
            assert(s>=0 && s<=1, 's must be in [0,1] range')
            assert(size(p,2)==3, 'p must be n x 3')
            n = length(p);
            i = fix(n*s)+1;
            dt = 1/(n-1);
            assert(i<=n-2,'cannot conpute 2nd derivative on last 2 points, take smaller s')
            p_s = p(i,:)';
            PP = ( p(i+1,:)' - p(i,:)' ) / dt ;
            PPP = ( p(i+2,:)' - 2 * p(i+1,:)' + p(i,:)' ) / dt ^ 2 ;
            t = PP / norm(PP);
            n = PPP / norm(PPP);
            b = cross(t,n);
            k = norm(cross(PP,PPP))/norm(PP)^3;
        end
        
        function draw_frenet(p,s)
            % help function to draw frenet frame
            [t,n,b,k,p_s]= SR_ROBO.frenet(p,s);
            clf;
            hold on
            title('FRENET FRAME AT S = '+ string(s) +', CURVATURE K = '+ string(k))
            plot3(p(:,1),p(:,2),p(:,3), 'r')
            plot3(p_s(1),p_s(2),p_s(3), 'k*')
            text(p_s(1),p_s(2),p_s(3), 'p(s)')
            t_s=p_s+t;
            n_s=p_s+n;
            b_s=p_s+b;
            quiver3(p_s(1),p_s(2),p_s(3),t(1),t(2),t(3),'r','LineWidth',2)
            quiver3(p_s(1),p_s(2),p_s(3),n(1),n(2),n(3),'g','LineWidth',2)
            quiver3(p_s(1),p_s(2),p_s(3),b(1),b(2),b(3),'b','LineWidth',2)
            text(t_s(1),t_s(2),t_s(3), 't(s)','Color',[1, 0 ,0])
            text(n_s(1),n_s(2),n_s(3), 'n(s)','Color',[0, 1 ,0])
            text(b_s(1),b_s(2),b_s(3), 'b(s)','Color',[0, 0 ,1])
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            view([110 15]);
            grid on
            pbaspect([1 1 1]);
            hold off
        end
        
        function [R_theta_s,theta_s,r,theta_ab] = axis_angle_trajectory(R_a,R_b,s,varargin)
            % R_theta_s; axis direction at varying of s
            % theta_s: varying of theta
            % r: axis
            % theta_ab: target angle
            assert(size(s,1)==1 && size(s,2)==1, 'theta(s) must be scalar')
            assert(s>=0 && s<=1, 's must be in [0,1] range')
            assert(size(R_a,1)==3 && size(R_a,2)==3 && size(R_b,1)==3 && size(R_b,2)==3, 'rotation matrices must be 3 x 3')
            [theta_ab,r]= SR_ROBO.ax_ang_inv(R_a'*R_b);
            assert(length(varargin)<=9, 'one condition of {linear,bcb,cubic,quintic,cubic_spline} must be specified')
            assert(isa(r,'double'),'no solution for the given axis angle combination')
            if ~isempty(varargin)
                if strcmp('linear', varargin{1})
                    assert(length(varargin)==1, '{linear} not requires more parameters')
                    linear = 1;
                elseif strcmp('bcb', varargin{1})
                    assert(length(varargin)==5,'{bcb} requires {"v_max",#v_max,"a_max",#a_max} parameters')
                    assert(isequal(varargin{2},'v_max'), 'first parameter should be v_max')
                    v_max = varargin{3};
                    assert(isa(v_max,'double'), 'v_max value should be double')
                    assert(isequal(varargin{4},'a_max'), 'second parameter should be a_max')
                    a_max = varargin{5};
                    assert(isa(a_max,'double'), 'a_max value should be double')
                    bcb = 1;
                elseif strcmp('cubic', varargin{1})
                    assert(length(varargin)==5,'{bcb} requires {"q_dot_0",#q_dot_0,"q_dot_1",#q_dot_1} parameters')
                    assert(isequal(varargin{2},'q_dot_0'), 'first parameter should be q_dot_0')
                    q_dot_0 = varargin{3};
                    assert(isa(q_dot_0,'double'), 'q_dot_0 value should be double')
                    assert(isequal(varargin{4},'q_dot_1'), 'second parameter should be q_dot_1')
                    q_dot_1 = varargin{5};
                    assert(isa(q_dot_1,'double'), 'q_dot_1 value should be double')
                    cubic = 1;
                elseif strcmp('quintic', varargin{1})
                    assert(length(varargin)==9,'{bcb} requires {"q_dot_0",#q_dot_0,"q_dot_1",#q_dot_1,"q_ddot_0",#q_ddot_0,"q_ddot_1",#q_ddot_1} parameters')
                    assert(isequal(varargin{2},'q_dot_0'), 'first parameter should be q_dot_0')
                    q_dot_0 = varargin{3};
                    assert(isa(q_dot_0,'double'), 'q_dot_0 value should be double')
                    assert(isequal(varargin{4},'q_dot_1'), 'second parameter should be q_dot_1')
                    q_dot_1 = varargin{5};
                    assert(isa(q_dot_1,'double'), 'q_dot_1 value should be double')
                    assert(isequal(varargin{6},'q_ddot_0'), 'third parameter should be q_ddot_0')
                    q_ddot_0 = varargin{7};
                    assert(isa(q_ddot_0,'double'), 'q_ddot_0 value should be double')
                    assert(isequal(varargin{8},'q_ddot_1'), 'fourth parameter should be q_ddot_1')
                    q_ddot_1 = varargin{9};
                    assert(isa(q_ddot_1,'double'), 'q_ddot_1 value should be double')
                    quintic = 1;
                elseif strcmp('cubic_spline', varargin{1})
                    assert(length(varargin)==3,'{cubic_spline} requires {"knots",{knots}} parameters')
                    assert(isequal(varargin{2},'knots'), 'first parameter should be v_max')
                    knots = varargin{3};
                    cubic_spline = 1;
                end
            end
            
            if exist('linear','var')
                theta_s = SR_ROBO.ptp_linear(0,theta_ab,s);
            elseif exist('bcb','var')
                [~,T,~] = SR_ROBO.ptp_bcb(theta_ab,a_max,v_max,s);
                if s == 0
                    fprintf('Real T = %f',T);
                end
                s = s * T;
                [theta_s,~,~] = SR_ROBO.ptp_bcb(theta_ab,a_max,v_max,s);
            elseif exist('cubic','var')
                theta_s = SR_ROBO.ptp_cubic(0,theta_ab,q_dot_0,q_dot_1,s);
            elseif exist('quintic','var')
                theta_s = SR_ROBO.ptp_quintic(0,theta_ab,q_dot_0,q_dot_1,q_ddot_0,q_ddot_1,s);
            elseif exist('cubic_spline','var')
                knots{1}{2}=0;
                knots{length(knots)}{2}=theta_ab;
                coeff = SR_ROBO.cubic_spline_coefficients(knots);
                theta_s = SR_ROBO.cubic_spline(coeff,knots,s);
            else
                theta_s = theta_ab * s;
            end
            R_theta_s = R_a * SR_ROBO.ax_ang(theta_s, r);
        end
        function draw_axis_angle_trajectory(R_theta_s_list)
            % help function to draw rotating reference frames 
            assert(size(R_theta_s_list,1)==3 && size(R_theta_s_list,2)==3, 'R_theta_s_list should be 3 x 3 x n')
            for i=1:size(R_theta_s_list,3)
                clf;
                title('AXIS ANGLE ROTATION TRAJECTORY')
                hold on
                quiver3(0,0,0,R_theta_s_list(1,1,i),R_theta_s_list(2,1,i),R_theta_s_list(3,1,i),'r','LineWidth',2)
                quiver3(0,0,0,R_theta_s_list(1,2,i),R_theta_s_list(2,2,i),R_theta_s_list(3,2,i),'g','LineWidth',2)
                quiver3(0,0,0,R_theta_s_list(1,3,i),R_theta_s_list(2,3,i),R_theta_s_list(3,3,i),'b','LineWidth',2)
                legend('x', 'y', 'z');
                pbaspect([1 1 1]);
                view([110 15]);
                xlim([-1 1]);
                ylim([-1 1]);
                zlim([-1 1]);
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                grid on
                drawnow
                hold off
            end
        end
    end
end
