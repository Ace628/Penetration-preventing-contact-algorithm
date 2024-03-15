function y=ufi_bar_EA_GA_PairForce_sym(x,rot,R0,E2,G,I2,A2)

Q = E2/G+A2*R0^2/I2+1;

C4 = -(E2*A2*E2*I2*R0+E2*I2*G*A2*R0+E2*A2*G*A2*R0^3)/(4*E2*A2*E2*I2*G*A2);
C5 = R0^3/pi/E2/I2;

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;
end

y =-2*(E2*A2+G*A2-G*A2*Q)/(G*A2*R0*Q)*C4*cos(X) - C5*(X - pi/2)/R0;

if rot > 0
    if x>=-pi+rot && x<= rot
        y = -y;
    end
elseif rot < 0
    if x>=-pi && x<=rot
        y = -y;
    elseif x >= pi+rot && x<=pi
        y = -y;
    end
elseif rot == 0
    if x>=-pi && x<=0
        y = -y;
    end
end

end
