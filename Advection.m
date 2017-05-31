close all;
%initial conditions for rho
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

%initial conditions for numerical approximation

t=0; tfinal=pi;
dx=x(2)-x(1); dy=y(2)-y(1); dt=0.2*dx;
rho=rho_0;
rho_x=zeros(N,M);
rho_y=zeros(N,M);

while t<tfinal
    if t+dt>tfinal
        dt=tfinal-t;
    end
    
    for i=1:N %looping over the grid
        for j=1:M
            %set the velocity for particular (i,j,t)
            c1=-cos(x(i))*sin(y(j))*cos(t);
            c2=sin(x(i))*cos(y(j))*cos(t);
            %boundary conditions
            if( i==1 || i==N || j==1 || j==M )
                rho(i,j)=0;
            else
                    
                if c1>0 
                    rho_x(i,j)=(rho_0(i,j)-rho_0(i-1,j));
                        if c2>0
                             rho_y(i,j)=(rho_0(i,j)-rho_0(i,j-1));

                        else
                             rho_y(i,j)=(rho_0(i,j+1)-rho_0(i,j));

                        end

                else 
                    rho_x(i,j)=(rho_0(i+1,j)-rho_0(i,j));
                        if c2>0
                             rho_y(i,j)=(rho_0(i,j)-rho_0(i,j-1));

                        else
                             rho_y(i,j)=(rho_0(i,j+1)-rho_0(i,j));

                        end
                 end
            end
            
            %THIS IS THE APPROXIMATION
            rho(i,j)=rho_0(i,j)-c1*dt/dx*rho_x(i,j)-c2*dt/dy*rho_y(i,j);
                       
        end
    end
    
    %prepare for next time step
    rho_0=rho; t=t+dt;
    %plot rho    
    mesh(x,y,rho);
    axis([-2 2 -2 2 0 1.1]);
    pause(dt);
end