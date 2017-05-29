close all;
%initial conditions
N=100;M=100;
center=[1 0];
r=0.25;
rho_0=zeros(N,M);
x=linspace(-pi/2,pi/2,N);
y=linspace(-pi/2,pi/2,M);

for i=1:N
    for j=1:M
        if (x(i)-center(1))^2+(y(j)-center(2))^2 < r^2
            rho_0(i,j)=1;
        end
    end
end

%other initial stuff
c1=1;c2=1;
t=0; tfinal=5;
dx=x(2)-x(1); dy=y(2)-y(1); dt=0.5*dx/abs(c1);
rho=rho_0;

%code

while t<tfinal
    if t+dt>tfinal
        dt=tfinal-t;
    end
    
    for i=2:N
        for j=2:M
        rho(i,j)=rho_0(i,j)-c1*dt/dx*(rho_0(i,j)-rho_0(i-1,j))-c2*dt/dy*(rho_0(i,j)-rho_0(i,j-1));
        %periodic boundary condition
            rho(1,j)=rho_0(1,j)-c1*dt/dx*(rho_0(1,j)-rho_0(N-1,j))-c2*dt/dy*(rho_0(1,j)-rho_0(N,j-1));
            rho(i,1)=rho_0(i,1)-c1*dt/dx*(rho_0(i,1)-rho_0(i-1,M))-c2*dt/dy*(rho_0(i,1)-rho_0(i,M-1));
            rho(1,1)=rho_0(1,1)-c1*dt/dx*(rho_0(1,1)-rho_0(N-1,M))-c2*dt/dy*(rho_0(1,1)-rho_0(N,M-1));
        end
    end
    
      
    %prepre for next time step
    rho_0=rho; t=t+dt;
        
    mesh(x,y,rho);
    axis([-2 2 -2 2 0 1.1]);
    pause(5*dt);
end
