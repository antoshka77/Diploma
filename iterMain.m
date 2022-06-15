%max levels
num_levels = 300;

%channel_count
channel_count = 1;

%Iaeneiaeuiia eiee?anoai ninoiyiee a y?ona
max_states_in_level = channel_count + 1;

%Ia?aue y?on n iaeneiaeuiui eiee?anoaii ninoiyiee
level_w_max_states = channel_count + 1;

%three first moments for entry and exit distribution
eps = 0.00001;
q_1_entr=1.4284837;
q_2_entr=3.3096436;
q_3_entr=8.9371421;
q_1_ex=     1;
q_2_ex=0.66666667;
q_3_ex=0.37037037;

%find Cox params
cox_req = CoxDistr;
if abs(q_1_entr * q_1_entr - q_2_entr) <= 10^(-4)
  'prostoy'
  cox_req.l1 = 1 / q_1_entr;
  cox_req.l2 = cox_req.l1;
  cox_req.y = 0.0;
else
  sum_first = q_3_entr - q_1_entr * q_2_entr;
  sum_second = 4 * q_1_entr^3 * q_3_entr + 4 * q_2_entr^3 + q_3_entr^2 - 3 * q_1_entr^2 * q_2_entr^2 - 6 * q_1_entr * q_2_entr * q_3_entr;
  denom = 2 * (q_2_entr - q_1_entr ^ 2);
  x_1 = (sum_first + sum_second^0.5)/denom;
  x_2 = (sum_first - sum_second^0.5)/denom;
  cox_req.l1 = 1/x_1;
  1/x_1
  cox_req.l2 = 1/x_2;
  1/x_2
  cox_req.y = (q_1_entr - x_1) / x_2;
  (q_1_entr - x_1) / x_2
end

fid = fopen('Cox_params.txt', 'w');
fprintf(fid, 'Entry:\n', q_1_entr);
fprintf(fid, 'cox_req.l1=%6.4g;\n', cox_req.l1);
fprintf(fid, 'cox_req.l2=%6.4g;\n', cox_req.l2);
fprintf(fid, 'cox_req.y=%6.4g;\n', cox_req.y);


cox_serv = CoxDistr;
if abs(q_1_ex * q_1_ex - q_2_ex) <= 10^(-4)
  'prostoy'
  cox_serv.l1 = 1 / q_1_ex;
  cox_serv.l2 = cox_serv.l1;
  cox_serv.y = 0.0;
else
  sum_first = q_3_ex - q_1_ex * q_2_ex;
  sum_second = 4 * q_1_ex^3 * q_3_ex + 4 * q_2_ex^3 + q_3_ex^2 - 3 * q_1_ex^2 * q_2_ex^2 - 6 * q_1_ex * q_2_ex * q_3_ex;
  denom = 2 * (q_2_ex - q_1_ex ^ 2);
  x_1 = (sum_first + sum_second^0.5)/denom;
  x_2 = (sum_first - sum_second^0.5)/denom;
  cox_serv.l1 = 1/x_1;
  1/x_1
  cox_serv.l2 = 1/x_2;
  1/x_2
  cox_serv.y = (q_1_ex - x_1) / x_2;
  (q_1_ex - x_1) / x_2
end

fprintf(fid, 'Exit:\n', q_1_ex);
fprintf(fid, 'cox_serv.l1=%6.4g;\n', cox_serv.l1);
fprintf(fid, 'cox_serv.l2=%6.4g;\n', cox_serv.l2);
fprintf(fid, 'cox_serv.y=%6.4g;\n', cox_serv.y);
fclose(fid);

%start approximation
start_approx_x = ((((1 - cox_req.y) * cox_req.l1 + cox_req.y * cox_req.l2) ...
    / ((1 - cox_serv.y) * cox_serv.l1 + cox_serv.y * cox_serv.l2)) ^ 4)^(1/5);

%matrix A
A = cell(num_levels,1);
for q = 1 : level_w_max_states - 1
    t = q + 1;
    help_arr = zeros(2 * q, 2 * t);
    for j = 1 : q
        help_arr(j, j) = (1 - cox_req.y) * cox_req.l1;
        help_arr(j + q, j) = cox_req.l2;
    end
    A{q} = help_arr;
end
help_arr = zeros(2 * level_w_max_states, 2 * level_w_max_states);
for q = 1 : level_w_max_states
    help_arr(q, q) = (1 - cox_req.y) * cox_req.l1;
    help_arr(q + level_w_max_states, q) = cox_req.l2;
end
for q = level_w_max_states : num_levels
    A{q} = help_arr;
end

%matrix B
B = cell(num_levels,1);
for q = 1 : level_w_max_states - 1
    help_arr = zeros(2 * (q + 1), 2 * q);
    for j = 1 : q
        diag = (q - (j - 1)) * cox_serv.l1 * (1 - cox_serv.y);
        help_arr(j, j) = diag;
        help_arr(j + 1, j) = j * cox_serv.l2;
        help_arr(q + j + 1, j + q) = diag;
        help_arr(q + j + 2, j + q) = j * cox_serv.l2;
    end
    B{q} = help_arr;
end
help_arr = zeros(2 * level_w_max_states, 2 * level_w_max_states);
q = level_w_max_states - 1;
for j = 1 : level_w_max_states
    diag = (q - (j - 1)) * cox_serv.l1 * (1 - cox_serv.y);
    help_arr(j, j) = diag;
    help_arr(j + level_w_max_states, j + level_w_max_states) = diag;
    if j ~= level_w_max_states
      help_arr(j + 1, j) = j * cox_serv.l2;
      help_arr(j + level_w_max_states + 1, j + level_w_max_states) = j * cox_serv.l2;
    endif
end
for q = level_w_max_states : num_levels
    B{q} = help_arr;
end

%matrix C
C = cell(num_levels,1);
for q = 1 : level_w_max_states
    help_arr = zeros(2 * q, 2 * q);
    for j = 1 : q - 1
        help_arr(j, j + 1) = (q - j) * cox_serv.l1 * cox_serv.y;
        help_arr(j, j + q) = cox_req.y * cox_req.l1;
        help_arr(j + q, j + q + 1) = help_arr(j, j + 1);
    end
    help_arr(q, 2 * (q)) = cox_req.y * cox_req.l1;
    C{q} = help_arr;
end
for q = level_w_max_states + 1 : num_levels
    C{q} = C{level_w_max_states};
end

%matrix D
D = cell(num_levels,1);
for q = 1: level_w_max_states
    help_arr = zeros(2 * q, 2 * q);
    sum_A = sum(A{q}, 2);
    sum_C = sum(C{q}, 2);
    if q > 1
        sum_B = sum(B{q - 1}, 2);
    else
        sum_B = [0; 0];
    end
    for j = 1 : 2 * q
        help_arr(j, j) = sum_A(j) + sum_B(j) + sum_C(j);
    end
    D{q} = help_arr;
end
for q = level_w_max_states + 1 : num_levels
    D{q} = D{level_w_max_states};
end

max_iterations=100;


%vector's of probablity initialization
t = cell(max_iterations, 1);
v = cell(1, num_levels + 1);
for q = 1 : level_w_max_states
    help_arr = zeros(1, 2 * q);
    for j = 1 : q
        help_arr(j) = 1 / (2 * q);
        help_arr(j + q) = 1 / (2 * q);
    end
    t{1}{q} = help_arr;
end
for q = level_w_max_states + 1 : num_levels - 1
    t{1}{q} = t{1}{level_w_max_states};
end
for q = 1 : level_w_max_states
    v{q} = ones(2 * q, 1);
end
for q = level_w_max_states + 1 : num_levels + 1
    v{q} = v{level_w_max_states};
end


beta1 = cell(max_iterations, 1);
beta2 = cell(max_iterations, 1);
c = cell(max_iterations, 1);
x = cell(max_iterations, 1);
z = cell(max_iterations, 1);
prob_final = cell(max_iterations, 1);

p = zeros(1, num_levels);

fin = 0;

%first approximation
Q = ((1/start_approx_x)*A{num_levels} + start_approx_x*B{num_levels})*inv(D{num_levels}-C{num_levels});
E = eye(2 * max_states_in_level);
Q = Q - E;
Q(1,:) = ones(1, 2 * max_states_in_level);
delta = zeros(2 * max_states_in_level, 1);
delta(1) = 1;
t{1}{num_levels} = (inv(Q)*delta)';

for q = 1: num_levels
    x{1}{q} = 0;
end;

fid = fopen('results_t.txt', 'w');
fid_2 = fopen('results_x.txt', 'w');
fid_3 = fopen('results_p_final.txt', 'w');
fid_4 = fopen('results_average_indicators.txt', 'w');
fprintf(fid, 'T values \nIteration number \n');
fprintf(fid_2, 'X values \nIteration number       Level %d ----> levels with smaller numbers ----> level 0 \n', num_levels);
fprintf(fid_3, 'Final probabilities of system states for %d levels and precision = %6.4g\n', num_levels, eps);
fprintf(fid_4, 'Final avarage indicators\n');
last_i = 0;
for i = 2:max_iterations

    % 1st iteration, last level
    beta1{i}{num_levels} = t{i-1}{num_levels - 1}*A{(num_levels - 1)}*inv(D{num_levels}-C{num_levels});
    beta2{i}{num_levels} = t{i-1}{num_levels}*B{num_levels}*inv(D{num_levels}-C{num_levels}); % t{1}{5} ili t{1}{4}, kak v knige? t{i-1}{num_levels-1}
    c{i}{num_levels} = (t{i-1}{num_levels-1}*B{num_levels}*v{(num_levels - 1)} - beta2{i}{num_levels}*A{num_levels}*v{num_levels})/(beta1{i}{num_levels}*A{num_levels}*v{(num_levels + 1)}); % same here  x{i}{num_levels} = 1/((c{i}{num_levels}*beta1{i}{num_levels} + beta2{i}{num_levels})*v{num_levels});
    x{i}{num_levels} = 1/((c{i}{num_levels}*beta1{i}{num_levels} + beta2{i}{num_levels})*v{num_levels});
    fprintf(fid_2, '%d       %6.4f ', (i-1), x{i}{num_levels}); %Print iteration number
    z{i}{num_levels} = c{i}{num_levels}*x{i}{num_levels};
    t{i}{num_levels} = z{i}{num_levels}*beta1{i}{num_levels} + x{i}{num_levels}*beta2{i}{num_levels};
    fprintf(fid, '%d    t{%d) = ', (i-1), num_levels); %Print iteration number
    fprintf(fid, '%6.4f ', t{i}{num_levels});
    fin = 0;

    %I?iaa?ea oneiaey ieii?aiey eoa?aoee
    if (( x{i}{num_levels} - x{i-1}{num_levels} ) < eps ) && ( ( x{i}{num_levels} - x{i-1}{num_levels} ) > 0 )
        fin = 1;
    end;

     % 1st iteration, levels all except num_levels and first
    for q = 1 : (num_levels - 2)

        j = num_levels - q;
        beta1{i}{j} = t{i-1}{j-1}*A{j-1}*inv(D{j}-C{j});
        beta2{i}{j} = t{i}{j+1}*B{j}*inv(D{j}-C{j});
        c{i}{j} = (t{i}{j+1}*B{j}*v{j} - beta2{i}{j}*A{j}*v{j+1})/(beta1{i}{j}*A{j}*v{j+1});
        x{i}{j} = 1/((c{i}{j}*beta1{i}{j} + beta2{i}{j})*v{j});
        fprintf(fid_2, '%6.4f ', x{i}{j});
        z{i}{j} = c{i}{j}*x{i}{j};
        t{i}{j} = z{i}{j}*beta1{i}{j} + x{i}{j}*beta2{i}{j};
        fprintf(fid, 't{%d) = ', j);
        fprintf(fid, '%6.4f ', t{i}{j});

        % I?iaa?ea oneiaey ieii?aiey eoa?aoee
        if (fin == 1) && (( x{i}{num_levels} - x{i-1}{num_levels} ) < eps ) && ( ( x{i}{num_levels} - x{i-1}{num_levels} ) > 0 )
            fin = 1;
        else
            fin = 0;
        end;
    end;

    %1st iteration, first level;

    beta2{i}{1} = t{i}{2}*B{1}*inv(D{1}-C{1});
    x{i}{1} = 1/(beta2{i}{1}*v{1});
    t{i}{1} = x{i}{1}*beta2{i}{1};
    fprintf(fid_2, '%6.4f  \n', x{i}{1});
    fprintf(fid, 't{1) = ');
    fprintf(fid, '%6.4f  \n', t{i}{1});


    if fin == 1
        last_i = i;
        break;
    end;

end;
if last_i == 0
  last_i = max_iterations;
end

%find vectors of final probablities
tmp = num_levels;
tmp_2 = 1;

for j = 1 : (num_levels - 1)
    for w = 1:j
        tmp_2 = tmp_2*x{last_i}{w};
    end;
    tmp = tmp + (num_levels - j - 1)*tmp_2;
    tmp_2 = 1;
end;
a = 1 / cox_req.l1 + cox_req.y / cox_req.l2;
b = 1 / cox_serv.l1 + cox_serv.y / cox_serv.l2;
p(1) = (num_levels - b/a)/tmp;
for j = 1 : (num_levels - 1)
    p(j+1) = p(j)*x{last_i}{j};
end;


fprintf(fid_2, '\nAbsolute probabilities of system being at a level\n');
total_p = 0;
for q = 0 : (num_levels - 1)
    j = num_levels - q;
    fprintf(fid_2, 'p(%d) = %6.4g ', j, p(j));
    total_p = total_p + p(j);
    prob_final{j} = t{last_i}{j}*p(j);
end;
fprintf(fid_2, ' \nprobability total = %6.4f ', total_p);


fprintf(fid_3, 'Level \n');
for i = 1:(level_w_max_states - 1)
    fprintf(fid_3, '%d ', i);
    for j = 1:2*i
        fprintf(fid_3, '%8.4g  ', prob_final{i}(j));
    end;
    fprintf(fid_3, ' \n ');
end;

for i = level_w_max_states : num_levels
    fprintf(fid_3, '%d ', i);
    for j = 1 : 2 * max_states_in_level
        fprintf(fid_3, '%8.4g  ', prob_final{i}(j));
    end;
    fprintf(fid_3, ' \n ');
end;


P_eff_ul1 = 0;
P_eff_l2 = 0;


for i = 1:(level_w_max_states - 1)
    fprintf(fid_3, '%d ', i);
    for j = 1:i
        P_eff_ul1 = P_eff_ul1 + prob_final{i}(j);
        if j ~= 1
          P_eff_l2 = P_eff_l2 + prob_final{i}(j + i);
        end;
    end;
end;

for i = level_w_max_states : num_levels - 1
    for j = 1 : max_states_in_level
      P_eff_ul1 = P_eff_ul1 + prob_final{i}(j);
      P_eff_l2 = P_eff_l2 + prob_final{i}(j + max_states_in_level);
    end;
end;

if P_eff_ul1 > 1
  P_eff_ul1 = 1.0;
end;
if P_eff_l2 > 1
  P_eff_l2 = 1.0;
end;
if P_eff_ul1 + P_eff_l2 > 1
  P_eff_l2 = 1 - P_eff_ul1;
end;

l_eff_ul1 = P_eff_ul1 * (1 - cox_req.y) * cox_req.l1;
l_eff_l2 = P_eff_l2 * cox_req.l2;
l_eff = l_eff_ul1 + l_eff_l2;


%averages calculating
n_average = 0;
n_average_q = 0;

for i = 1 : num_levels
    n_average = n_average + p(i) * (i - 1);
    if i > level_w_max_states
      n_average_q = n_average_q + p(i) * (i - level_w_max_states);
    endif
end;

intens_average = 1 / q_1_entr;
fprintf(fid_4, 'intens_average = %6.4f\n', intens_average);
fprintf(fid_4, 'n_average = %6.4f\n', n_average);
fprintf(fid_4, 't_average = %6.4f\n', n_average / intens_average);
fprintf(fid_4, 'n_average_queue = %6.4f\n', n_average_q);
fprintf(fid_4, 't_average_queue = %6.4f\n', n_average_q / intens_average); %???
fprintf(fid_4, 'n_average_ch = %6.4f\n', n_average - n_average_q); %???




%Cae?uoea oaeeia
fclose(fid);
fclose(fid_2);
fclose(fid_3);
fclose(fid_4);

%Vitt's formula
[n_average, n_average_q, n_average_c] = get_averages(1/q_1_entr, q_2_entr*2, 1/q_1_ex, 2*q_2_ex, max_states_in_level-1)
t_average = n_average / intens_average
t_average_q = n_average_q / intens_average



