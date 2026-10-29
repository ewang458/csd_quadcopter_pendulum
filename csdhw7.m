s = tf('s');
P = 10 / (s*(s + 10));

A_border = 2;
B_border = 300;

[mag_dist, ~] = bode(P, A_border);
[mag_noise, ~] = bode(P, B_border);


k_min = (1/0.02 - 1) / squeeze(mag_dist);
k_max = 0.05 / squeeze(mag_noise);


k = 180;
L1 = k * P;

[Gm, Pm, Wcg, Wcp] = margin(L1);

fprintf('Gain crossover frequency : %.2f \n', Wcp);
fprintf('Phase Margin: %.2f \n', Pm);

figure;
margin(L1);
grid on;