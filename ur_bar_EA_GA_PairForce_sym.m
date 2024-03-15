function y=ur_bar_EA_GA_PairForce_sym(x,rot,R0,E2,G,I2,A2)

C1 = R0/2/G/A2;
C4 = -(E2*A2*E2*I2*R0+E2*I2*G*A2*R0+E2*A2*G*A2*R0^3)/(4*E2*A2*E2*I2*G*A2);
C5 = R0^3/pi/E2/I2;

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;
end

y = C1*sin(X) +C4*(sin(X)-(X - pi/2)*cos(X)) + C5;

end
