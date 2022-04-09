P_tan = date_indiv_SS(225);

H = [P_tan.den{1}(2) P_tan.den{1}(4) 0; P_tan.den{1}(1) P_tan.den{1}(3) 0; 0 P_tan.den{1}(2) P_tan.den{1}(4)];
det1 = det(H(1, 1));
det2 = det(H(1:2, 1:2));
det3 = det(H);
numitor = P_tan.den{1};
poli = roots(numitor);

t = (0:0.01:180)';
h_pondere = impulse(P_tan,t);
rasp_trp = step(P_tan,t);
trp= double(t>=0);
rasp_conv = conv(trp,h_pondere)*0.01;
rasp_conv=rasp_conv(1:size(t));
norm_dif=norm(rasp_trp - rasp_conv,inf);

x0=[1 1 1];
rasp_tot=lsim(ss_ci(P_tan),trp,t,x0);
rasp_perm=evalfr(P_tan*trp,0);
rasp_tran=rasp_tot-rasp_perm;
rasp_libr=initial(ss_ci(P_tan),x0,t);
rasp_fort=rasp_tot-rasp_libr;

tc1=stepinfo(P_tan).RiseTime;
tt1=stepinfo(P_tan).SettlingTime;
tv1=stepinfo(P_tan).PeakTime;
sr1=stepinfo(P_tan).Overshoot;
s=tf('s');
P_aux= 1/(10*s+1);
tc2=stepinfo(P_tan*P_aux).RiseTime;
tt2=stepinfo(P_tan*P_aux).SettlingTime;
tv2=stepinfo(P_tan*P_aux).PeakTime;
sr2=stepinfo(P_tan*P_aux).Overshoot;
tc3=stepinfo(P_tan*(tf('s')+1)).RiseTime;
tt3=stepinfo(P_tan*(tf('s')+1)).SettlingTime;
tv3=stepinfo(P_tan*(tf('s')+1)).PeakTime;
sr3=stepinfo(P_tan*(tf('s')+1)).Overshoot;

save('tema_225.mat', 'H', 'det1', 'det2', ...
'det3', 'poli', 'h_pondere', 'rasp_trp', 'rasp_conv', ...
'norm_dif', 'rasp_tot', 'rasp_perm', 'rasp_tran', ...
'rasp_libr', 'rasp_fort', 'tc1', 'tt1', 'tv1', 'sr1', ...
'tc2', 'tt2', 'tv2', 'sr2', 'tc3', 'tt3', 'tv3', 'sr3');
