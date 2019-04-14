function r=analytical(x)
  pi = 3.14159;
  r = (-2 * sqrt(2*pi) .* (-1 + 8 * pi.^2 * x.^2)) ./ exp(2 * pi^2 .* x.^2);
endfunction
