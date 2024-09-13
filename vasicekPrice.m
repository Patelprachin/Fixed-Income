function price = vasicekPrice(alpha,r,sigma,R0,maturity,decimal,timeNow,coupon)
    % Note that
    C = @(T,t) (1-exp(-alpha*(T-t)))/alpha; % part of the pricing function
    A = @(T,t) r*(T-t)-r*(1-exp(-alpha*(T-t)))/alpha-sigma^2*(1-exp(-2*alpha*(T-t)))/...
        (4*alpha^3)+sigma^2*(1-exp(-alpha*(T-t)))/(alpha^3)-sigma^2*(T-t)/(2*alpha^2); % part of the pricing function
    
    if maturity > 1
        % if the maturity > 1, the payoff will be 1+coupon, else 1.
        q1 = coupon;
    else
        q1 = 0;
    end
    
    t = timeNow; % addresses timeNow-value, this is always set to zero
    
    if decimal > 0
        % price for the first payoff.
        price = q1 * exp(-(R0*C(decimal,t) + A(decimal,t)));
    else
        price = 0;
    end
    
    for i = t + 1:maturity - 1
        price = price + q1 * exp(-(R0*C(i + decimal,t) + A(i + decimal,t)));
    end
    
    price = price + (1 + q1) * exp(-(R0*C(maturity + decimal,t) + A(maturity + decimal,t)));
end
