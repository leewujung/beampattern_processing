function E = get_ellipse_param(coef)

% Equation coef
a = coef(1);
b = coef(2);
c = coef(3);
d = coef(4);
f = coef(5);
g = coef(6);
E.coef = coef;
E.eqt = sprintf('%2.5f*x^2 + %2.5f*x*y + %2.5f*y^2 + %2.5f*x+ %2.5f*y + %2.5f',a,b,c,d,f,g);

% Semi-axes
E.a0 = sqrt((2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g))/((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
E.b0 = sqrt((2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g))/((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));
E.ar = E.a0/E.b0;
if E.ar<1
    E.ar = 1/E.ar;
end

% Counterclockwise angle from +x-axis to the major axis
if b==0
    if a<c
        E.theta = 0;
    else
        E.theta = pi/2;
    end
else
    if a<c
        E.theta = 1/2*acot((a-c)/(2*b));
    else
        E.theta = pi/2+1/2*acot((a-c)/(2*b));
    end
end

% Eccentricity (0<=e<1, e=0: circle)
E.e = sqrt(1-b^2/a^2);


