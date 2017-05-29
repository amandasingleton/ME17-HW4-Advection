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

t=0; tfinal=5;
dx=x(2)-x(1); dy=y(2)-y(1); dt=0.5*dx;
rho=rho_0;

%code

while t<tfinal
    if t+dt>tfinal
        dt=tfinal-t;
    end
    
    for i=2:N-1
        for j=2:M-1
            c1=-cos(i)*sin(j)*cos(dt);
            c2=sin(i)*cos(j)*cos(dt);
            if c1>0 %%this is probably going to break when zero
                rho_x(i,j)=(rho_0(i,j)-rho_0(i-1,j));
                rho(1,j)=rho_0(1,j)-c1*dt/dx*(rho_0(1,j)-rho_0(N-1,j))-c2*dt/dy*(rho_0(1,j)-rho_0(N,j-1));
            elseif c1<0 
                rho_x(i,j)=(rho_0(i+1,j)-rho_0(i,j));
                rho(1,j)=rho_0(1,j)-c1*dt/dx*(rho_0(2,j)-rho_0(N,j))-c2*dt/dy*(rho_0(1,j+1)-rho_0(N,j));
            end
            if c2>0
                rho_y(i,j)=(rho_0(i,j)-rho_0(i,j-1));
                rho(i,1)=rho_0(i,1)-c1*dt/dx*(rho_0(i,1)-rho_0(i-1,M))-c2*dt/dy*(rho_0(i,1)-rho_0(i,M-1));
            elseif c2<0
                rho_y(i,j)=(rho_0(i,j+1)-rho_0(i,j));
                rho(i,1)=rho_0(i,1)-c1*dt/dx*(rho_0(i+1,1)-rho_0(i,M))-c2*dt/dy*(rho_0(i,2)-rho_0(i,M));
            end
            rho(i,j)=rho_0(i,j)-c1*dt/dx*rho_x(i,j)-c2*dt/dy*rho_y(i,j);
            %periodic boundary condition
            
            
            rho(1,1)=rho_0(1,1)-c1*dt/dx*(rho_0(1,1)-rho_0(N-1,M))-c2*dt/dy*(rho_0(1,1)-rho_0(N,M-1));
        end
    end
    
      
    %prepre for next time step
    rho_0=rho; t=t+dt;
        
    mesh(x,y,rho);
    %axis([-2 2 -2 2 0 1.1]);
    %view(0, 90);
    pause(5*dt);
end
