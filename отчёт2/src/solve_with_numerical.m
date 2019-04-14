function F=solve_with_numerical(f, a, b, M)
  F = zeros(M, 1);
  us = linspace(-b, b, M);
  
  for u_index = 1:M
    u = us(u_index);
    F(u_index) = num_integration(@(x) f(x) .* exp(-2 .* pi .* i .* u .* x), -a, a, 512);
  endfor
endfunction
