function r=extend_with_zeros(arr, number_of_zeros)
  n = length(arr);
  r = zeros(2 * number_of_zeros + n, 1);
  
  for i=1:n
    r(i + number_of_zeros) = arr(i);
  endfor
endfunction
