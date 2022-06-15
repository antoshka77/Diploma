function [n_average, n_average_q, n_average_c] = get_averages(l_entr, l_2_entr, l_ex, l_2_ex, chanel_count)
  ro = l_entr / (l_ex * chanel_count)
  coef = 1;
  for i=1:chanel_count
    coef = coef * (ro * chanel_count);
  end
  ch_factor = factorial(chanel_count);
  sum_1 = 1;
  for i=1:chanel_count-1
    coef_1 = 1;
    for j=1:i
      coef_1 = coef_1 * (ro * chanel_count);
    endfor
    sum_1 = sum_1 + coef_1 / factorial(i);
  end
  P = (coef / (ch_factor * (1 - ro))) * (1 / (sum_1 + coef / (ch_factor * (1 - ro))));
  n_average_q = ro * P / (1 - ro);
  n_average = ro * chanel_count * (1 + P / (chanel_count * (1 - ro)));
  n_average_c = ro * chanel_count;
  D_entr = l_2_entr - 1 / (l_entr * l_entr);
  D_ex = l_2_ex - 1 / (l_ex * l_ex);
  c_entr = D_entr * (l_entr * l_entr);
  c_ex = D_ex * (l_ex * l_ex);
  coef = (c_entr + c_ex) / 2;
  n_average_q = n_average_q * coef;
  n_average_c = n_average_c * coef;
  n_average = n_average * coef;
end

