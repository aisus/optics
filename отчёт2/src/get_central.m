function r=get_central(arr, num_zeros, n)
  r = zeros(n, 1);
  for i = 1:n
    r(i) = arr(i + num_zeros);
  endfor
endfunction
