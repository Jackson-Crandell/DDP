function [dynamics] = fnDynamics(m_c, m_p, l, g)
x = sym('x', [1 4]); % State vector.
u = sym('u');

den = (m_c + m_p * sin(x(3)) ^ 2);
Fa = [
    x(2);
    (m_p * sin(x(3)) * (l * (x(4) ^ 2) + g * cos(x(3)))) / den;
    x(4);
    (- m_p * l * (x(4) ^ 2) * cos(x(3)) * sin(x(3)) - (mc + m_p) * g * sin(x(3))) / (l * den)
];

Fb = [0; u / den; 0; (-u * cos(x(3))) / (l * den)];

% System dynamics F(x,u).
F = Fa + Fb;

dynamics = struct();
dynamics.F = matlabFunction(F);
dynamics.Fx = matlabFunction(jacobian(F, x));
dynamics.Fu = matlabFunction(jacobian(F, u));
dynamics.Fb = matlabFunction(Fb); % Used to com_pute added noise due to stochastic controls.
end