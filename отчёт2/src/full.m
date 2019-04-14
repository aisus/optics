
function r=analytical(x)
  pi = 3.14159;
  r = (-2 * sqrt(2*pi) .* (-1 + 8 * pi.^2 * x.^2)) ./ exp(2 * pi^2 .* x.^2);
endfunction

function r=extend_with_zeros(arr, number_of_zeros)
  n = length(arr);
  r = zeros(2 * number_of_zeros + n, 1);
  
  for i=1:n
    r(i + number_of_zeros) = arr(i);
  endfor
endfunction

function r=gauss_f(x)
  r = exp(-x .^ 2);
endfunction

function r=get_central(arr, num_zeros, n)
  r = zeros(n, 1);
  for i = 1:n
    r(i) = arr(i + num_zeros);
  endfor
endfunction

function r=my_f(x)
  r = (4 .* (x .^ 2) .- 2) .* exp(-x .^ 2 ./ 2);
endfunction

function integral = num_integration(f,a,b,n)
    width = (b-a)/n;
    x = linspace(a,b,n);
    integral = width * sum( f( (x(1:n-1)+x(2:n))/2 ) );
end

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

function F=solve_with_numerical(f, a, b, M)
  F = zeros(M, 1);
  us = linspace(-b, b, M);
  
  for u_index = 1:M
    u = us(u_index);
    F(u_index) = num_integration(@(x) f(x) .* exp(-2 .* pi .* i .* u .* x), -a, a, 512);
  endfor
endfunction

function r=swap_parts(arr)
  n = length(arr);
  half = floor(n / 2);
  r = vertcat(arr(half + 1:n), arr(1:half));
endfunction
