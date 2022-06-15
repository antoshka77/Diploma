func_1 = @(x, l) l*exp(-x*l); %M
func_2 = @(x, l, k) (l^k)*exp(-x*l).*(x.^(k-1))/factorial(k-1); %Ek
func_3 = @(x, a_1, l_1, a_2, l_2) (a_1 * func_1(x, l_1) + a_2 * func_1(x, l_2)); %H2

l_M_0_7 = 0.7;
l_M_1_0 = 1.0;
l_M_2_1 = 2.1;

l_E3_0_7 = 2.1;
l_E3_1_0 = 3.0;
l_E3_2_1 = 6.3;


l_1_H2_0_7 = 0.3561;
l_2_H2_0_7 = 1.9659;
a_1_H2_0_7 = 0.4;
a_2_H2_0_7 = 0.6;
l_1_H2_1_0 = 2/3;
l_2_H2_1_0 = 1.5;
a_1_H2_1_0 = 0.4;
a_2_H2_1_0 = 0.6;
l_1_H2_2_1 = 1.0683;
l_2_H2_2_1 = 5.8707;
a_1_H2_2_1 = 0.4;
a_2_H2_2_1 = 0.6;

q_1_D_0_7 = 10/7;
q_2_D_0_7 = q_1_D_0_7 ^ 2;
q_3_D_0_7 = q_1_D_0_7 ^ 3;
q_1_D_1_0 = 1.0;
q_2_D_1_0 = q_1_D_1_0 ^ 2;
q_3_D_1_0 = q_1_D_1_0 ^ 3;
q_1_D_2_1 = 1/2.1;
q_2_D_2_1 = q_1_D_2_1 ^ 2;
q_3_D_2_1 = q_1_D_2_1 ^ 3;

fid = fopen('integrate_moments.txt', 'w');

q_1=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x, 0, Inf);
q_2=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x.^2, 0, Inf);
q_3=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

k=3;
q_1=integral(@(x)func_2(x, l_E3_1_0, k).*x, 0, Inf);
q_2=integral(@(x)func_2(x, l_E3_1_0, k).*x.^2, 0, Inf);
q_3=integral(@(x)func_2(x, l_E3_1_0, k).*x.^3, 0, Inf);
fprintf(fid, 'q_1_ex=%6.8g;\n', q_1);
fprintf(fid, 'q_2_ex=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_ex=%6.8g;\n', q_3 / 6);

fclose(fid);


%M entry ch=1
q_1=integral(@(x)func_1(x, l_M_0_7).*x, 0, Inf);
q_2=integral(@(x)func_1(x, l_M_0_7).*x.^2, 0, Inf);
q_3=integral(@(x)func_1(x, l_M_0_7).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%M entry ch=3
q_1=integral(@(x)func_1(x, l_M_2_1).*x, 0, Inf);
q_2=integral(@(x)func_1(x, l_M_2_1).*x.^2, 0, Inf);
q_3=integral(@(x)func_1(x, l_M_2_1).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%M service
q_1=integral(@(x)func_1(x, l_M_1_0).*x, 0, Inf);
q_2=integral(@(x)func_1(x, l_M_1_0).*x.^2, 0, Inf);
q_3=integral(@(x)func_1(x, l_M_1_0).*x.^3, 0, Inf);
fprintf(fid, 'q_1_ex=%6.8g;\n', q_1);
fprintf(fid, 'q_2_ex=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_ex=%6.8g;\n', q_3 / 6);

%E3 entry ch=1
k=3;
q_1=integral(@(x)func_2(x, l_E3_0_7, k).*x, 0, Inf);
q_2=integral(@(x)func_2(x, l_E3_0_7, k).*x.^2, 0, Inf);
q_3=integral(@(x)func_2(x, l_E3_0_7, k).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%E3 entry ch=3
k=3;
q_1=integral(@(x)func_2(x, l_E3_2_1, k).*x, 0, Inf);
q_2=integral(@(x)func_2(x, l_E3_2_1, k).*x.^2, 0, Inf);
q_3=integral(@(x)func_2(x, l_E3_2_1, k).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%E3 service
k=3;
q_1=integral(@(x)func_2(x, l_E3_1_0, k).*x, 0, Inf);
q_2=integral(@(x)func_2(x, l_E3_1_0, k).*x.^2, 0, Inf);
q_3=integral(@(x)func_2(x, l_E3_1_0, k).*x.^3, 0, Inf);
fprintf(fid, 'q_1_ex=%6.8g;\n', q_1);
fprintf(fid, 'q_2_ex=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_ex=%6.8g;\n', q_3 / 6);

%H2 entry ch=1
q_1=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x, 0, Inf);
q_2=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x.^2, 0, Inf);
q_3=integral(@(x)func_3(x, a_1_H2_0_7, l_1_H2_0_7, a_2_H2_0_7, l_2_H2_0_7).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%H2 entry ch=3
q_1=integral(@(x)func_3(x, a_1_H2_2_1, l_1_H2_2_1, a_2_H2_2_1, l_2_H2_2_1).*x, 0, Inf);
q_2=integral(@(x)func_3(x, a_1_H2_2_1, l_1_H2_2_1, a_2_H2_2_1, l_2_H2_2_1).*x.^2, 0, Inf);
q_3=integral(@(x)func_3(x, a_1_H2_2_1, l_1_H2_2_1, a_2_H2_2_1, l_2_H2_2_1).*x.^3, 0, Inf);
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3 / 6);

%H2 service
q_1=integral(@(x)func_3(x, a_1_H2_1_0, l_1_H2_1_0, a_2_H2_1_0, l_2_H2_1_0).*x, 0, Inf);
q_2=integral(@(x)func_3(x, a_1_H2_1_0, l_1_H2_1_0, a_2_H2_1_0, l_2_H2_1_0).*x.^2, 0, Inf);
q_3=integral(@(x)func_3(x, a_1_H2_1_0, l_1_H2_1_0, a_2_H2_1_0, l_2_H2_1_0).*x.^3, 0, Inf);
fprintf(fid, 'q_1_ex=%6.8g;\n', q_1);
fprintf(fid, 'q_2_ex=%6.8g;\n', q_2 / 2);
fprintf(fid, 'q_3_ex=%6.8g;\n', q_3 / 6);

%D entry ch=1
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1_D_0_7);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2_D_0_7 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3_D_0_7 / 6);

%D entry ch=3
fprintf(fid, 'q_1_entr=%6.8g;\n', q_1_D_2_1);
fprintf(fid, 'q_2_entr=%6.8g;\n', q_2_D_2_1 / 2);
fprintf(fid, 'q_3_entr=%6.8g;\n', q_3_D_2_1 / 6);

%D service
fprintf(fid, 'q_1_ex=%6.8g;\n', q_1_D_1_0);
fprintf(fid, 'q_2_ex=%6.8g;\n', q_2_D_1_0 / 2);
fprintf(fid, 'q_3_ex=%6.8g;\n', q_3_D_1_0 / 6);
