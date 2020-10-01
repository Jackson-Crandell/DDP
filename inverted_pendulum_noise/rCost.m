function  [l0,l_x,l_xx,l_u,l_uu,l_ux] = rCost(x, u, k,R,P,target, dt)

l0 = 0.5*(u'*R*u + (x - target)'*P*(x - target));
l_x = P*(x - target);
l_xx = P; 

l_u = R * u;
l_uu = R;
l_ux = zeros(1,2);


% l0 = u'*R*u;
% l_x = zeros(2,1);
% l_xx = zeros(2,2); 
% l_u = R * u;
% l_uu = R;
% l_ux = zeros(1,2);
