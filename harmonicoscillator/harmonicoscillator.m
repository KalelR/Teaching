
w = 2;
z0 = [5,5];
nu = 0;
t0 = 0;
tfinal = 15;
t = t0:0.01:tfinal;
Z1 = z0(1) .* cos(w .* t) + (z0(2)/w) .* sin(w.*t);
Z2 = z0(2) .* cos(w .* t) - (w*z0(1)) .* sin(w.*t);
Z = [Z1; Z2];

figure()
plot(t, Z)

figure()
plot(Z(1,:), Z(2,:))


y0 = z0;  
p = [nu,w];
[t,y] = ode23(@(t,y)oscillator(t,y,p),[t0 tfinal], y0);
figure()
plot(t,y)

function yp = oscillator(t,y, p)
    nu = p(1); omega = p(2);
    yp = [y(2), -2*nu*y(2) - omega^2 * y(1)]';
end

