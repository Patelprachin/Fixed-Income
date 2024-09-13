
function [ri,rall] = VasicekModelCW(a , mu, sigma,r0,T,nsteps,M)
    dt = T/nsteps;
    Vts = (sigma^2*(1-exp(-2*a*dt)))/(2*a);
    rall = [];
    for j = 1:M
	r = zeros(nsteps+1,1); r(1) = r0;
	Z = randn(1,nsteps);
	for i = 1:nsteps
        Its = r(i)*exp(-a*dt) + mu*(1-exp(-a*dt));
        r(i+1) = Its + sqrt(Vts)*Z(i);
   end
	rall = [rall, r];
    end
    ri = [rall(:,end)] ;
end