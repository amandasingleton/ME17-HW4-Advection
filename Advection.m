close all;
%initial conditions for rho
N=100;M=100;
center=[1 0];
r=0.25;
rho_0=zeros(N,M);
rho_x=zeros(N,M);
rho_y=zeros(N,M);
x=linspace(-pi/2,pi/2,N);
y=linspace(-pi/2,pi/2,M);

cx = @(x1,y1,t1) -cos(x1)*sin(y1)*cos(t1); 
cy = @(x1,y1,t1) sin(x1)*cos(y1)*cos(t1);
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

while t<tfinal
    if t+dt>tfinal
        dt=tfinal-t;
    end
    
    for i=2:N-1
        for j=2:M-1
            c1=cx(x(i),y(j),t);
            c2=cy(x(i),y(j),t);
            if( i==1 || i==N || j==1 || j==M )
                rho(i,j)=0;
            end
            if c1>0 
                rho_x(i,j)=(rho_0(i,j)-rho_0(i-1,j));
                    if c2>0
                         rho_y(i,j)=(rho_0(i,j)-rho_0(i,j-1));
                         
                    else
                         rho_y(i,j)=(rho_0(i,j+1)-rho_0(i,j));
                         
                    end
                                    
            elseif c1<0 
                rho_x(i,j)=(rho_0(i+1,j)-rho_0(i,j));
                    if c2>0
                         rho_y(i,j)=(rho_0(i,j)-rho_0(i,j-1));
                        
                    else
                         rho_y(i,j)=(rho_0(i,j+1)-rho_0(i,j));
                         
                    end
            end
           
            rho(i,j)=rho_0(i,j)-c1*dt/dx*rho_x(i,j)-c2*dt/dy*rho_y(i,j);
              
                       
        end
    end
    
    %prepre for next time step
    rho_0=rho; t=t+dt;
        
    mesh(x,y,rho);
    axis([-2 2 -2 2 0 1.1]);
    xlabel('x');
    %view(0, 90);
    pause(dt);
end