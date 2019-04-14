function integral = num_integration(f,a,b,n)
    width = (b-a)/n;
    x = linspace(a,b,n);
    integral = width * sum( f( (x(1:n-1)+x(2:n))/2 ) );
end
