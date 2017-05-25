N=100;M=100;
center=[1 0];
r=0.25;
rho=zeros(N,M);
x=linspace(-pi/2,pi/2,N);
y=linspace(-pi/2,pi/2,M);

for i=1:N
    for j=1:M
        if (x(i)-center(1))^2+(y(j)-center(2))^2 < r^2
            rho(i,j)=1;
        end
    end
end

mesh(x,y,rho);