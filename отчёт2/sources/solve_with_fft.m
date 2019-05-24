function [F, b]=solve_with_fft(f, a, N, M)
  xs = transpose(linspace(-a, a, N + 1))(1:N);
  h_x = 2 * a / N;

  discrete_f = f(xs);
  num_zeros = floor((M - N) / 2);
  extended_f = extend_with_zeros(discrete_f, num_zeros);
  swapped_f = swap_parts(extended_f);
  
  F = fft(swapped_f) * h_x;
  F = swap_parts(F);
  F = get_central(F, num_zeros, N);
  
  b = N^2 / (4 * a * M);
endfunction
