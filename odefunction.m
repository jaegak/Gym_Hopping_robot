function dX=odefunction(t,X)
    %% Set the variable
    global dx Torque
    x0=X(1);
    dx0=X(2);
    x1=X(3);
    dx1=X(4);
    x2=X(5);
    x2=0.75*t;
    dx2=X(6);
    y0=X(7);
    dy0=X(8);
    y1=X(9);
    dy1=X(10);
    y2=X(11);
    dy2=X(12);
    q1=X(13);
    dq1=X(14);
    q2=X(15);
    dq2=X(16);
    w=X(17);
    dw=X(18);
    x=X(19);
    P_gain = 5;
    %% limitation
    if w>1.2
        w=1.2;
        if dw>0
            dw=0;
        end
    end
    if w<0.4
        w=0.4;
        if dw<0
            dw=0;
        end
    end
    if x>w-0.4
        x=w-0.4;
        if dx>0
            dx=0;
        end
    end
    
    if x>0.12
        x=0.12;
        if dx>0
            dx=0;
        end
    end
    if x<0.05
        x=0.05;
        if dx<0
            dx=0;
        end
    end
    %% Touch down setting
    global xTD TD btTD ctTD by0 bdy0 cH Top k H bH
    if bH ~= H
        cH = 0;
        bH = H;
    end
    
    if TD==0
        if y0<=0
            xTD=x0;
            btTD=ctTD;
            ctTD=t;
            TD=1;
            Top=0;
            k=0;
        end
    end
    if TD==1
        if y0>0
            TD=0;
        end
    end
    
%     if Top==0
        if y0 > cH
%             if bdy0>0
%                 if dy0<0
                    cH = y0;
%                     k=k+1;
%                     if k == 1
%                         Top=1;
%                     end
%                 end
%             end
        end
%     end
    bdy0=dy0;    
    by0=y0;
    %% Appendix 2 Simulation Parameters
    M1 = 1 ; % kg
    M2 = 10 ; % kg
    J1 = 1 ; % kg-m^2
    J2 = 10 ; % kg-m^2
    r1 = 0.5 ; % m
    r2 = 0.4 ; % m
    k0 = 1 ; % m
    K_L = 10^3 ; % Nt/m
    K_L2 = 10^5 ; % Nt/m
    B_L2 = 125 ; % Nt-s/m
    K_G = 10^4 ; % Nt/m
    B_G = 75 ; %Nt-s/m
%     H = 0.4 ; % m
    K_P=1200 ;  % Nt-m/rad
    if y0<=0
        K_P=1800;
    end
    
    K_V=60; % Nt-m-s/rad
    if y0<=0
        K_V=200;
    end
    
    g=9.81; % m/s^2
    
    %% Appendix 1 Equation of Motion for Motel of Fig 1
    W=w-r1;

    %% eq45

    F_K=K_L2*(k0-w+x)-B_L2*dw;
    if (k0-w+x) > 0
        F_K=K_L*(k0-w+x);
    end
    
    %% eq46
    
    F_x=0;
    if y0<0        
        F_x = -K_G*(x0-xTD)-B_G*dx0; %-xTD
    end
    
    %% eq47
    F_y=0;
    if y0<0
        F_y = -K_G*y0-B_G*dy0;
    end
    
    A=zeros(5,5);
    B=zeros(5,1);
    
    %% eq 40
    A(1,1)=cos(q1)*(M2*W*w+J1);
    A(1,2)=M2*r2*W*cos(q2);
    A(1,3)=M2*W;
    A(1,4)=0;
    A(1,5)=M2*W*sin(q1);
    
%     B(1,1)=W*M2*(dq1^2*W*sin(q1)-2*dq1*dw*cos(q1)+...
%         r2*dq2^2*sin(q2)+r1*dq1^2*sin(q1))-...
%         r1*F_x*cos(q1)^2+cos(q1)*(r1*F_y*sin(q1)-Torque)+...
%         F_K*W*sin(q1);
    B(1,1)=M2*W*(w*dq1^2*sin(q1) + r2*dq2^2*sin(q2) - w*dw*dq1*cos(q1)) -...
        F_x*r1*cos(q1)^2 + F_y*r1*cos(q1)*sin(q1) - Torque*cos(q1) + F_K*W*sin(q1);
    %% eq 41
    A(2,1)=-sin(q1)*(M2*W*w+J1);
    A(2,2)=-M2*r2*W*sin(q2);
    A(2,3)=0;
    A(2,4)=M2*W;
    A(2,5)=M2*W*cos(q1);
    
%     B(2,1)=W*M2*(dq1^2*W*cos(q1)+2*dq1*dw*sin(q1)+r2*dq2^2*cos(q2)+...
%         r1*dq1^2*cos(q1)-g)+...
%         r1*F_x*cos(q1)*sin(q1)-sin(q1)*...
%         (r1*F_y*sin(q1)-Torque)+F_K*W*cos(q1);
    B(2,1)=M2*W*(w*dq1^2*cos(q1) + r2*dq2^2*cos(q2) + 2*dw*dq1*sin(q1) - g) +...
        F_x*r1*cos(q1)*sin(q1) - F_y*r1*sin(q1)^2 + Torque*sin(q1) + F_K*W*cos(q1);
    %% eq 42
    A(3,1)=cos(q1)*(M1*r1*W-J1);
    A(3,2)=0;
    A(3,3)=M1*W;
    A(3,4)=0;
    A(3,5)=0;
    
%     B(3,1)=W*(M1*r1*dq1^2*sin(q1) - F_K*sin(q1)+F_x)-...
%         cos(q1)*(F_y*r1*sin(q1)-F_x*r1*cos(q1)-Torque);
    B(3,1)= W*(F_x - F_K*sin(q1) + M1*r1*dq1^2*sin(q1)) - cos(q1)*(-F_x*r1*cos(q1) + F_y*r1*sin(q1) - Torque);
    %% eq 43
    A(4,1)=-sin(q1)*(M1*r1*W-J1);
    A(4,2)=0;
    A(4,3)=0;
    A(4,4)=M1*W;
    A(4,5)=0;
    
%     B(4,1)=W*(M1*r1*dq1^2*cos(q1) - F_K*cos(q1)+F_y-M1*g)+...
%         sin(q1)*(r1*F_y*sin(q1)-r1*F_x*cos(q1)-Torque);
    B(4,1)= W*(F_y - F_K*cos(q1) - M1*g + M1*r1*dq1^2*cos(q1)) + sin(q1)*(-F_x*r1*cos(q1) + F_y*r1*sin(q1) - Torque);
    %% eq 44
    A(5,1)=-cos(q2-q1)*J1*r2;
    A(5,2)=J2*W;
    A(5,3)=0;
    A(5,4)=0;
    A(5,5)=0;
    
    B(5,1)=W*(F_K*r2*sin(q2-q1)+Torque)-r2*cos(q2-q1)*(r1*F_y*sin(q1)-r1*F_x*cos(q1)-Torque);
    
    %% eq 40~44
    C=A\B;
    
    ddq1=C(1);
    ddq2=C(2);
    ddx0=C(3);
    ddy0=C(4);
    ddw=C(5);
        
    %% eq 30~33
    ddy1=ddy0-r1*(ddq1*sin(q1)+dq1^2*cos(q1));
    ddx1=ddx0+r1*(ddq1*cos(q1)-dq1^2*sin(q1));
    
    ddy2=ddy0+ddw*cos(q1)-w*ddq1*sin(q1)-w*dq1^2*cos(q1)-r2*(ddq2*sin(q2)+dq2^2*cos(q2))-2*dw*dq1*sin(q1);
    ddx2=ddx0+ddw*sin(q1)+w*ddq1*cos(q1)-w*dq1^2*sin(q1)+r2*(ddq2*cos(q2)-dq2^2*sin(q2))+2*dw*dq1*cos(q1);
    
    %% Set the result    
    dX(1)=dx0;
    dX(2)=ddx0;
    dX(3)=dx1;
    dX(4)=ddx1;
    dX(5)=dx2;
    dX(6)=ddx2;
    dX(7)=dy0;
    dX(8)=ddy0;
    dX(9)=dy1;
    dX(10)=ddy1;
    dX(11)=dy2;
    dX(12)=ddy2;
    dX(13)=dq1;
    dX(14)=ddq1;
    dX(15)=dq2;
    dX(16)=ddq2;
    dX(17)=dw;
    dX(18)=ddw;
    dX(19)=dx*P_gain;
    dX=dX';
end
