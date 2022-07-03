% created by Mohsen Lotfi 
% Persian Gulf University
%-------------------------------------------------------------------------%

% Simulating the 1-D viscous Burgers' equation 
% Numerical scheme used is a FOU & SOU & QUICK
%-------------------------------------------------------------------------%

%%
%Specifying Parameters for Courant-Num = 0.7 & t=2 & 10s & 20s & 30S &40s
clear all
clc

% FOR MESH WITH 20 CELL

nx=20;                    %Number of steps in space(x)
dt=0.06842;               %Width of each time step 
nt=146;                   %Number of time steps for  t=10
% nt=292;                   %Number of time steps for  t=20
% nt=438;                   %Number of time steps for  t=30
% nt=585;                   %Number of time steps for  t=40

% FOR MESH WITH 40 CELL:
 
% nx=40;                    %Number of steps in space(x)
% dt=0.033333;              %Width of each time step 
% nt=300;                   %Number of time steps for  t=10
% nt=600;                   %Number of time steps for  t=20
% nt=900;                   %Number of time steps for  t=30
% nt=1200;                  %Number of time steps for  t=40

% FOR MESH WITH 80 CELL:

nx=80;                    %Number of steps in space(x)
dt=0.0165;                %Width of each time step 
% nt=606;                   %Number of time steps for  t=10
nt=1212;                  %Number of time steps for  t=20
% nt=1818;                  %Number of time steps for  t=30
% nt=2424;                  %Number of time steps for  t=40

% FOR MESH WITH 160 CELL:

% nx=160;                   %Number of steps in space(x)
% dt=0.0082;                %Width of each time step
% nt=1219;                  %Number of time steps for  t=10
% nt=2439;                  %Number of time steps for  t=20
% nt=3658;                  %Number of time steps for  t=30
% nt=4878;                  %Number of time steps for  t=40

dx=20/(nx-1);               %Width of space step
x=0:dx:20;                  %Range of x (0,2) and specifying the grid points
u=zeros(1,nx);              %Preallocating u
un=zeros(1,nx);             %Preallocating un
CN=abs(dt/dx);              %Courant Number

%%
%Initial conditions

for i=1:nx-1
    if (0<=x(i) && x(i)<=1)
       u(i)=0;
    elseif (1<=x(i) && x(i)<=3)
               u(i)=0.5*(x(i)-1);
        elseif (3<=x(i) && x(i)<=5)
               u(i)=0.5*(5-x(i));
    elseif (5<=x(i) && x(i)<=6)
                    u(i)=0;
    end
    x(i+1)=x(i)+dx;
end

% ui=u;

%%
% Boundary condition
u(1)=0;
u(nx)=0;

%%
% plotting Numerical solution:
i=3:nx-1;
for it=1:nt
    un=u;  
    
% FOU (FTBS Explicit):
   
%         u(i)=un(i)-(CN/2)*((un(i)).^2-(un(i-1)).^2);

% SOU:

%         u(i)=un(i)-(CN/2)*((3/2*un(i)-1/2*un(i-1)).^2-((3/2)*un(i-1)-(1/2)*un(i-2)).^2);
        
% QUICK: 

%         u(i)=un(i)- (CN/2)*(((3/8)*un(i+1)+(3/4)*un(i)-(1/8)*un(i-1)).^2-((3/8)*un(i)+(3/4)*un(i-1)-(1/8)*un(i-2)).^2) ;

        
% LAX_WENDROFF

        u(i)=un(i)-(CN/2)*((1/2)*(un(i+1)).^2-(1/2)*(un(i-1)).^2)+(12*(CN.^2)).*((1/2)*(un(i)+un(i+1)).*((1/2)*((un(i+1)).^2-(1/2).*(un(i).^2)))-((1/2)*(un(i)+un(i-1))).*((1/2)*((un(i)).^2-(1/2)*((un(i-1)).^2))));
          
% LAX_FRIEDRIHS

%         u(i)=(1/2)*(un(i-1)+un(i+1))-(CN/4)*(((un(i+1)).^2-(un(i-1)).^2));
end    

%     figure(3)
%     h=plot(x,u,'--*b');       %plotting the velocity profile
%     legend('t = 10')

    figure(4)
    h=plot(x,u,'--*r');       %plotting the velocity profile
    legend('t = 20')

%     figure(5)
%     h=plot(x,u,'--*g');       %plotting the velocity profile
%     legend('t = 30')

%     figure(6)
%     h=plot(x,u,'--*m');       %plotting the velocity profile
%     legend('t = 40')
 
% axis([0 20 0 1.1])
% title({'1-D Inviscous Burger Equation';['time(\itt) = ',num2str(dt*it)]})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('Transport property profile (u) \rightarrow')

save u

%%
% plotting the Exact solution:

t=dt*nt;
x0=2;
x=0;  
dx=20/nx; 
u=zeros(1,nx);
L=2*sqrt(x0*t/2);
H=sqrt(2*x0/t);

if t<x0
     
    for i=1:nx-1
    if x(i) < x0+t
    u(i)=x(i)/(x0+t);
    elseif x(i) >= x0+t
        u(i) =  1+ ((x0+t-x(i))/(x0-t));
    end

x(i+1) =x(i)+dx;
end

elseif t>=x0
    
    for i=1:nx-1
        if x(i)<= L
            u(i)= (x(i)/L)*H;
        elseif x(i) > L
            u(i)=0;

        end
    x(i+1) =x(i) +dx;
    end
    
end %if t<x0

u_exact=u;
u_exact2=0*u;
u_exact2((1/dx)+1:nx)=u_exact(1:nx-(1/dx));

% figure(1)
% plot(x,u_exact2)
% hold on
% axis([0 20 0 1.1])
% title({'Exact Burger Equation'})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('Transport property profile (u) \rightarrow')


% Calculating Error:
 
 load u
 Error = abs(u_exact2 - u);
 Err = (norm(Error)/norm(u_exact2))
