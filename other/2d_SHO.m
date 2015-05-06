counter = 1;
hbar = 1.06e-34;
m = 9.1e-31;
a = 1e-10;
q = 1.6e-19;
b = hbar^2/(2*m*a^2*q);
omega = 0.96e15;
Nx = 31;
Ny = Nx;

H = zeros(Nx*Ny, Nx*Ny);

for ix=1:Nx;
    for iy=1:Ny;
        NN= ix + (iy-1)*Nx;
        V= ( (a*(ix-1-(Nx-1)/2))^2 + (a*(iy-1-(Ny-1)/2))^2)...
              *(1/2)*m*(omega)^2/q;
        H(NN,NN)=4*b+V;
        if ix<Nx
            NN2=ix+1+(iy-1)*Nx;
            H(NN,NN2)=-b;
            H(NN2,NN)=-b;
        end
        if iy<Ny
            NN2=ix+iy*Nx;
            H(NN,NN2)=-b;
            H(NN2,NN)=-b;
        end
    end
end

[psis, evals] = eig(H);
psi1 = reshape(psis(:,1),Nx,Ny);
surface(psi1)

for E=0:0.1:1.4;
    k=acos(1-E/(2*b));
    Hme=H-E*eye(NN);
    Tau1=zeros(NN,1);
    Tau1(1,1)=-b;
    Tau2=zeros(NN,1);
    Tau2(NN,1)=-b;
    A1=zeros(NN+2);
    A1(1:NN,1:NN)=Hme;
    A1(1:NN,NN+1)=Tau1;
    A1(1:NN,NN+2)=Tau2;
    A1(NN+1,1:NN)=Tau1;
    A1(NN+2,1:NN)=Tau2;
    A1(NN+1,NN+1)=2*b-b*exp(i*k)-E;
    A1(NN+2,NN+2)=2*b-b*exp(i*k)-E;
    A2=zeros(NN+2,1);
    A2(1:NN,1)=-Tau1;
    A2(NN+1,1)=-2*b+b*exp(-i*k)+E;
    d=linsolve(A1,A2);
    psi(counter,:)=d(1:NN);
    R(counter)=abs(d(NN+1))^2;
    T(counter)=abs(d(NN+2))^2;
    counter=counter+1;
    eg=eig(H)/(hbar*omega/q);
end

plot(eg)