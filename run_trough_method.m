ch_n=1;
max_state_level_num = ch_n+1;
max_levels=300;
f=cell(1, max_levels+1);
F=cell(max_levels, 1);
F{1} = B{1}*inv(D{1}-C{1});
for i = 1:ch_n
  j = i+1;
  F{j} = B{j} * inv(D{j} - C{j} - F{j-1} * A{j-1});
end
j=ch_n+1;
for k=j:max_levels
  i = k+1;
  F{i} = B{j} * inv(D{j} - C{j} - F{i-1} * A{j});
end
Q = ((1/start_approx_x)*A{max_levels} + start_approx_x*B{max_levels})*inv(D{max_levels}-C{max_levels});
E = eye(2 * max_state_level_num);
Q = Q - E;
%Q(1,:) = ones(1, 2 * max_state_level_num+1);
delta = zeros(2 * max_state_level_num, 1);
delta(1) = 1;
f{num_levels+1} = (inv(Q)*delta)';
for k = max_levels:-1:1
  f{k} = f{k+1} * F{k};
end
sum_t_i = 0;
for i=1:max_levels+1
  sum_t_i = sum_t_i + sum(f{i});
end
p=cell(1, max_levels+1);
for i=1:max_levels+1
  f{i} = f{i} / sum_t_i;
  p{i} = sum(f{i});
end
n_average = 0;
n_average_q = 0;
for i = 1 : max_levels+1
    n_average = n_average + p{i} * (i - 1);
    if i > max_state_level_num
      n_average_q = n_average_q + p{i} * (i - max_state_level_num);
    endif
end;
n_average
t_average = n_average / intens_average
n_average_q
t_average_q = n_average_q / intens_average
n_average_ch = n_average - n_average_q
sum(cell2mat(p))
