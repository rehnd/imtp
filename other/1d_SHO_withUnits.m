% Same as main, but with units

hbar = 1.06e-34;
m = 9.1e-31;
a = 1e-12;
q = 1.6e-19;
b = hbar^2/(2*m*a^2*q);
omega = 0.96e15;
n = 1001;
N = n;

xrange = [-800*a 800*a];
h = (xrange(2) - xrange(1))/(N-1);
x = linspace(-N/2+1/2,N/2-1/2,N);

D = sparse(1:n,1:n, 2*b ,n,n);
E = sparse(2:n,1:n-1, -b*ones(1,n-1),n,n);
T = (E + D + E');

V = sparse(1:n,1:n, (a*x).^2*(1/2)*m*omega^2/q,n,n);

H = T + V;
H = full(H);
[psi, ev] = eig(H);

figure;
plot(a*x,psi(:,1:3))
hold on;
plot(a*x,diag(V))
ylim([-0.1,.1])
xlim(xrange)