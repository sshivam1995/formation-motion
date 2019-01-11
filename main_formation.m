close all
%function [xleaderl,xleader2,xfollower1,xfollower2,uleader,ufollower,rleader,rfollower]=reg(tf,dt,alpha)
%function reg(tf,dt,alpha)
tf=27; dt=0.01; alpha=30; alpha_tail=30;
% follower trajectory

speed_up=20;

a=0.2;
T=0.2;
deltat=0.01*T;

L=1.5;

A=[ 0 0 1 0;
    0 0 0 1;
    0 0 -a 0;
    0 0 0 -a];

B= [0 0;
    0 0;
    1 0;
    0 1];

sat=5;

% C= [0 0 1 0;
%     0 0 0 1];

C=[1 0 0 0;
   0 1 0 0];

t=0:dt:tf;

N=size(t,2);
Nbar=(N-1)*23/27;


xleader=zeros(4,N);
xfollower1=zeros(4,N);
xfollower2=zeros(4,N);
xfollower21=zeros(4,N);
xfollower22=zeros(4,N);
xfollower23=zeros(4,N);


X0_leader=[0;0;0;6*speed_up/tf];
X0_follower1=[-L/2;-sqrt(3)*L/2;0;6*speed_up/tf];
X0_follower2=[L/2;-sqrt(3)*L/2;0;6*speed_up/tf];
X0_follower21=[-2*L/2;-2*sqrt(3)*L/2;0;6*speed_up/tf];
X0_follower22=[0;-2*sqrt(3)*L/2;0;6*speed_up/tf];
X0_follower23=[2*L/2;-2*sqrt(3)*L/2;0;6*speed_up/tf];


xleader(:,1)=X0_leader;
xfollower1(:,1)=X0_follower1;
xfollower2(:,1)=X0_follower2;
xfollower21(:,1)=X0_follower21;
xfollower22(:,1)=X0_follower22;
xfollower23(:,1)=X0_follower23;

uleader=zeros(2,N);
ufollower1=zeros(2,N);
ufollower2=zeros(2,N);
ufollower21=zeros(2,N);
ufollower22=zeros(2,N);
ufollower23=zeros(2,N);

rvelleader=zeros(2,N);
rleader=zeros(2,N);
predictleader=zeros(2,N);

rvelfollower1=zeros(2,N);
rfollower1=zeros(2,N);
predictfollower1=zeros(2,N);

rvelfollower2=zeros(2,N);
rfollower2=zeros(2,N);
predictfollower2=zeros(2,N);

nleader=zeros(1,N);
nfollower1=zeros(1,N);
nfollower2=zeros(1,N);
nfollower21=zeros(1,N);
nfollower22=zeros(1,N);
nfollower23=zeros(1,N);


for i=1:(N-1)/27*10
    rvelleader(:,i)=[2*pi*sin(2*pi*i/(N*20/27))/tf,2*pi*cos(2*pi*i/(N*20/27))/tf]';
end

for i=(N-1)/27*10+1:(N-1)*20/27
    rvelleader(:,i)=[-2*pi*sin(2*pi*i/(N*20/27))/tf,2*pi*cos(2*pi*i/(N*20/27))/tf]';
end

for i=(N-1)*20/27+1:N
    rvelleader(:,i)=[0,2*pi/tf]';
end    

rvelleader=speed_up*rvelleader;

for i=2:N
    rleader(:,i)=rleader(:,i-1)+rvelleader(:,i-1)*dt;
end    

for i=1:N-T/dt
   predictleader(:,i)=rleader(:,i+T/dt); 
end

for i=N-T/dt+1:N
   predictleader(:,i)=predictleader(:,N-T/dt); 
end    

reffollower1(:,1)=C*(xleader(:,1)-xfollower1(:,1));
reffollower2(:,1)=C*(xleader(:,1)-xfollower2(:,1));

for i=2:N 
    xdot=dxdt(xleader(:,i-1),uleader(:,i-1),A,B);
    xleader(:,i)=xleader(:,i-1)+xdot*dt;
    
    [gu(:,i),guprime]=future_calc(xleader(:,i-1),uleader(:,i-1),A,B,C, T, deltat);
    
    deltau=guprime\(predictleader(:,i)-gu(:,i))*alpha*dt;
    uleader(:,i)=uleader(:,i-1)+deltau;
    
     % Begin  saturation u 
   % if norm(uleader(:,i))>1
   %     uleader(:,i)=1*uleader(:,i)/norm(uleader(:,i));
    % end
     %  end saturation
    
     
    
     xdot=dxdt(xfollower1(:,i-1),ufollower1(:,i-1),A,B);
     xfollower1(:,i)=xfollower1(:,i-1)+xdot*dt;     
     [gu1(:,i),guprime1]=future_calc(xfollower1(:,i-1),ufollower1(:,i-1),A,B,C, T, deltat);
     
     xdot=dxdt(xfollower2(:,i-1),ufollower2(:,i-1),A,B);
     xfollower2(:,i)=xfollower2(:,i-1)+xdot*dt;     
     [gu2(:,i),guprime2]=future_calc(xfollower2(:,i-1),ufollower2(:,i-1),A,B,C, T, deltat);
     
     
     xdot=dxdt(xfollower21(:,i-1),ufollower21(:,i-1),A,B);
     xfollower21(:,i)=xfollower21(:,i-1)+xdot*dt;     
     [gu21(:,i),guprime21]=future_calc(xfollower21(:,i-1),ufollower21(:,i-1),A,B,C, T, deltat);
     
     xdot=dxdt(xfollower22(:,i-1),ufollower22(:,i-1),A,B);
     xfollower22(:,i)=xfollower22(:,i-1)+xdot*dt;     
     [gu22(:,i),guprime22]=future_calc(xfollower22(:,i-1),ufollower22(:,i-1),A,B,C, T, deltat);
     
     xdot=dxdt(xfollower23(:,i-1),ufollower23(:,i-1),A,B);
     xfollower23(:,i)=xfollower23(:,i-1)+xdot*dt;     
     [gu23(:,i),guprime23]=future_calc(xfollower23(:,i-1),ufollower23(:,i-1),A,B,C, T, deltat);
     
    
    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xleader(1:2,i)-xleader(1:2,1);
%     else    
%         p=xleader(1:2,i)-xleader(1:2,i-ki+T/dt);
%     end    
        
        p=xleader(1:2,i)-xleader(1:2,i-1);
        

        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
%               q1=rot1*p;
%               q2=rot2*p;

        end
        
        reffollower1(:,i)=C*xleader(:,i)+L*q1+p*T/dt;
        pos21(:,i)=C*xleader(:,i)+L*q1;
        reffollower2(:,i)=C*xleader(:,i)+L*q2+p*T/dt;
        pos22(:,i)=C*xleader(:,i)+L*q2;
        
        pos31(:,i)=C*xleader(:,i)+2*L*q1;
        pos33(:,i)=C*xleader(:,i)+2*L*q2;
        pos32(:,i)=C*xleader(:,i)+L*q1+L*q2;
    
    deltau1=guprime1\(reffollower1(:,i)-gu1(:,i))*alpha*dt;
    ufollower1(:,i)=ufollower1(:,i-1)+deltau1;
    
    deltau2=guprime2\(reffollower2(:,i)-gu2(:,i))*alpha*dt;
    ufollower2(:,i)=ufollower2(:,i-1)+deltau2;
    
    
    
%%   for 3rd layer

%% follower 21 and 22_1
    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xfollower1(1:2,i)-xfollower1(1:2,1);
%     else    
%         p=xfollower1(1:2,i)-xfollower1(1:2,i-ki+T/dt);
%     end    
        
        p=xfollower1(1:2,i)-xfollower1(1:2,i-1);

        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
        end
        
        reffollower21(:,i)=C*xfollower1(:,i)+L*q1+p*T/dt;
        reffollower22_1(:,i)=C*xfollower1(:,i)+L*q2+p*T/dt;
  %      pos31(:,i)=C*xfollower1(:,i)+L*q1;
  %      pos32_1(:,i)=C*xfollower1(:,i)+L*q2;
        
%% follower 22_2 and 23

    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xfollower2(1:2,i)-xfollower2(1:2,1);
%     else    
%         p=xfollower2(1:2,i)-xfollower2(1:2,i-ki+T/dt);
%     end    

        
        p=xfollower2(1:2,i)-xfollower2(1:2,i-1);
        
        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
        end
        
        reffollower22_2(:,i)=C*xfollower2(:,i)+L*q1+p*T/dt;
        reffollower23(:,i)=C*xfollower2(:,i)+L*q2+p*T/dt;       
        
   %     pos32_2(:,i)=C*xfollower2(:,i)+L*q1;
    %    pos33(:,i)=C*xfollower2(:,i)+L*q2;
        
    %    pos32(:,i)=0.5*(pos32_2(:,i)+pos32_1(:,i));
        
    reffollower22(:,i)=(reffollower22_1(:,i)+reffollower22_2(:,i))/2;    
        
    deltau21=guprime21\(reffollower21(:,i)-gu21(:,i))*alpha_tail*dt;
    ufollower21(:,i)=ufollower21(:,i-1)+deltau21;
    
    deltau22=guprime22\(reffollower22(:,i)-gu22(:,i))*alpha_tail*dt;
    ufollower22(:,i)=ufollower22(:,i-1)+deltau22;
    
    deltau23=guprime23\(reffollower23(:,i)-gu23(:,i))*alpha_tail*dt;
    ufollower23(:,i)=ufollower23(:,i-1)+deltau23;
        
     % Begin  saturation u 
   % if norm(ufollower1(:,i))>3
   %     ufollower1(:,i)=3*ufollower1(:,i)/norm(ufollower1(:,i));
    % end
     %  end saturation
   
   
     % Begin  saturation u 
   % if norm(ufollower(:,i))>3
   %     ufollower(:,i)=3*ufollower(:,i)/norm(ufollower(:,i));
    % end
     %  end saturation
     
     %% add saturation
     
     if norm(uleader(:,i))>sat
       uleader(:,i)=sat*uleader(:,i)/norm(uleader(:,i));
     end
     
     if norm(ufollower1(:,i))>sat
       ufollower1(:,i)=sat*ufollower1(:,i)/norm(ufollower1(:,i));
     end
  
     if norm(ufollower2(:,i))>sat
       ufollower2(:,i)=sat*ufollower2(:,i)/norm(ufollower2(:,i));
     end
     
     if norm(ufollower21(:,i))>sat
       ufollower21(:,i)=sat*ufollower21(:,i)/norm(ufollower21(:,i));
     end
     
     if norm(ufollower22(:,i))>sat
       ufollower22(:,i)=sat*ufollower22(:,i)/norm(ufollower22(:,i));
     end
     
     if norm(ufollower23(:,i))>sat
       ufollower23(:,i)=sat*ufollower23(:,i)/norm(ufollower23(:,i));
     end
     
   
     nleader(i)=norm(uleader(:,i));
     nfollower1(i)=norm(ufollower1(:,i));
     nfollower2(i)=norm(ufollower2(:,i));
     nfollower21(i)=norm(ufollower21(:,i));
     nfollower22(i)=norm(ufollower22(:,i));
     nfollower23(i)=norm(ufollower23(:,i));
end



%% input norm graph
    figure(4)

    plot(t(1:Nbar),nleader(1:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(1:Nbar),nfollower1(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower1(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower21(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower22(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower23(1:Nbar), 'LineWidth',1.5);
    
    x1=xlabel('Time$~[s]$');
 y1=ylabel('Input norm$~[m/s^2]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('acceleration vs time')
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('Input_norm','-dsvg','-r0')


%% r-g

for i=1:N
r_g_leader(i)=norm(predictleader(:,i)-gu(:,i));
r_g_follower1(i)=norm(reffollower1(:,i)-gu1(:,i));
r_g_follower2(i)=norm(reffollower2(:,i)-gu2(:,i));
r_g_follower21(i)=norm(reffollower21(:,i)-gu21(:,i));
r_g_follower22(i)=norm(reffollower22(:,i)-gu22(:,i));
r_g_follower23(i)=norm(reffollower23(:,i)-gu23(:,i));
end

figure(42)

    plot(t(2:Nbar),r_g_leader(2:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(2:Nbar),r_g_follower1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower21(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower22(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower23(2:Nbar), 'LineWidth',1.5);


    
        x1=xlabel('Time$~[s]$');
 y1=ylabel('Control Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('Control error vs time')
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('r_g','-dsvg','-r0')

%% total error
for i=1:N
error_1_1(i)=norm(rleader(:,i)-C*xleader(:,i));
error_2_1(i)=norm(pos21(:,i)-C*xfollower1(:,i));
error_2_2(i)=norm(pos22(:,i)-C*xfollower2(:,i));
error_3_1(i)=norm(pos31(:,i)-C*xfollower21(:,i));
error_3_2(i)=norm(pos32(:,i)-C*xfollower22(:,i));
error_3_3(i)=norm(pos33(:,i)-C*xfollower23(:,i));
end

figure(21)

    plot(t(2:Nbar),error_1_1(2:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(2:Nbar),error_2_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_2_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_3(2:Nbar), 'LineWidth',1.5);


    
        x1=xlabel('Time$~[s]$');
 y1=ylabel('Total Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('Total Error vs time')
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('total_error','-dsvg','-r0')


%% inter_agent distance

figure(97)
for i=1:Nbar
   dL_11(i)= norm(xleader(1:2,i)-xfollower1(1:2,i));
   dL_12(i)= norm(xleader(1:2,i)-xfollower2(1:2,i));
   
   d11_21(i)= norm(xfollower1(1:2,i)-xfollower21(1:2,i));
   d11_22(i)= norm(xfollower1(1:2,i)-xfollower22(1:2,i));
   d12_22(i)= norm(xfollower2(1:2,i)-xfollower22(1:2,i));
   d12_23(i)= norm(xfollower2(1:2,i)-xfollower23(1:2,i));
end    

plot (t(1:Nbar),dL_11(1:Nbar), 'LineWidth',1.5);
hold on
plot (t(1:Nbar),dL_12(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d11_21(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d11_22(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d12_22(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d12_23(1:Nbar), 'LineWidth',1.5);

       x1=xlabel('Time$~[s]$');
 y1=ylabel('Inter-agent distance$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$d_{A_{1,1}-A_{2,1}}$','$d_{A_{1,1}-A_{2,2}}$','$d_{A_{2,1}-A_{3,1}}$','$d_{A_{2,1}-A_{3,2}}$','$d_{A_{2,2}-A_{3,2}}$','$d_{A_{2,2}-A_{3,3}}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('Inter-agent distance vs time')
hold off

ylim([1.35 1.65])

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('d_ij','-dsvg','-r0')

%% lateral error

figure (77)
for i=1:100
    r1(:,i)=[predictleader(1,1)-1+i/100;predictleader(2,1)];
    r2(:,i)=[reffollower1(1,1)-1+i/100;reffollower1(2,1)];
    r3(:,i)=[reffollower2(1,1)-1+i/100;reffollower2(2,1)];
    r4(:,i)=[reffollower21(1,1)-1+i/100;reffollower21(2,1)];
    r5(:,i)=[reffollower22(1,1)-1+i/100;reffollower22(2,1)];
    r6(:,i)=[reffollower23(1,1)-1+i/100;reffollower23(2,1)];
    
end    
    
for i=1:N
    r1(:,i+100)=predictleader(:,i);  
    r2(:,i+100)=reffollower1(:,i);
    r3(:,i+100)=reffollower2(:,i);
    r4(:,i+100)=reffollower21(:,i);
    r5(:,i+100)=reffollower22(:,i);
    r6(:,i+100)=reffollower23(:,i);
end



    
    [close_point,lateral_error_norm(1,:),arc]=distance2curve(r1',xleader(1:2,:)');
    [close_point,lateral_error_norm(2,:),arc]=distance2curve(r2',xfollower1(1:2,:)');
    [close_point,lateral_error_norm(3,:),arc]=distance2curve(r3',xfollower2(1:2,:)');
    [close_point,lateral_error_norm(4,:),arc]=distance2curve(r4',xfollower21(1:2,:)');
    [close_point,lateral_error_norm(5,:),arc]=distance2curve(r5',xfollower22(1:2,:)');
    [close_point,lateral_error_norm(6,:),arc]=distance2curve(r6',xfollower23(1:2,:)');
    
    

 plot(t(1:Nbar),lateral_error_norm(1,1:Nbar),'LineWidth',1.5)
 hold on 
 plot(t(1:Nbar),lateral_error_norm(2,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(3,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(4,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(5,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(6,1:Nbar),'LineWidth',1.5)

 
  x1=xlabel('Time$~[s]$');
 y1=ylabel('Lateral Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('Lateral error vs time')
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('lat_error','-dsvg','-r0')


 %% trajectories
    figure(13)
        
    plot(xleader(1,1:Nbar),xleader(2,1:Nbar),'LineWidth',1.5)
    hold on
    plot(xfollower1(1,1:Nbar),xfollower1(2,1:Nbar),'LineWidth',1.5)

    plot(xfollower2(1,1:Nbar),xfollower2(2,1:Nbar),'LineWidth',1.5)

    title('trajectories')
    plot(xfollower21(1,1:Nbar),xfollower21(2,1:Nbar),'LineWidth',1.5)

    plot(xfollower22(1,1:Nbar),xfollower22(2,1:Nbar),'LineWidth',1.5)

    plot(xfollower23(1,1:Nbar),xfollower23(2,1:Nbar),'LineWidth',1.5)
    
    pbaspect([2.5 1 1])
    
    
    xlim([-5 95]);
    ylim([-20 20]);
    
    
        x1=xlabel('$Z_1~[s]$');
 y1=ylabel('$Z_2~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
%     
% s1 = plot(xleader(1,1),xleader(2,1),'o','MarkerSize', 4,'MarkerFaceColor','red');      % bot 1
% %q = plot(r(1,1),r(2,1),'o','MarkerFaceColor','blue');                          % ref
% s2 = plot(xfollower1(1,1),xfollower1(2,1),'o','MarkerSize', 4,'MarkerFaceColor','green');
% s3 = plot(xfollower2(1,1),xfollower2(2,1),'o','MarkerSize', 4,'MarkerFaceColor','green');
% s4 = plot(xfollower21(1,1),xfollower21(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
% s5 = plot(xfollower22(1,1),xfollower22(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
% s6 = plot(xfollower23(1,1),xfollower23(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
%  
% for k = 2:Nbar
%     s1.XData = xleader(1,k);
%     s1.YData = xleader(2,k);
%       
% %     q.XData = r(1,k);
% %     q.YData = r(2,k);
%     
%     s2.XData = xfollower1(1,k);
%     s2.YData = xfollower1(2,k);
%     
%     s3.XData = xfollower2(1,k);
%     s3.YData = xfollower2(2,k);
% 
%     s4.XData = xfollower21(1,k);
%     s4.YData = xfollower21(2,k);
%     
%     
%     s5.XData = xfollower22(1,k);
%     s5.YData = xfollower22(2,k);
%     
%     
%     s6.XData = xfollower23(1,k);
%     s6.YData = xfollower23(2,k);
%     
%     
%     drawnow
% end


 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
  
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off


fig.PaperUnits = 'inches';
print('Path','-dsvg','-r0')










