vd = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
vg10_noqc_10nm=[0,3.59 * 1e-7, 4.14 * 1e-7,4.41 * 1e-7, 4.59 * 1e-7, 4.64 * 1e-7];
vg10_qc_10nm=[0,2.41 * 1e-7, 2.9 * 1e-7, 3.2 * 1e-7,3.3 * 1e-7, 3.4 * 1e-7];
vg06_noqc_10nm=[0,6.13 * 1e-8, 7.0 * 1e-8 ,7.7 * 1e-8,8.4 * 1e-8, 8.7 * 1e-8];
vg06_qc_10nm=[0,6.2 * 1e-8, 7.2 * 1e-8,8.0 * 1e-8, 8.6 * 1e-8, 9.0 * 1e-8];
vg08_noqc_10nm=[0,2.12 * 1e-7, 2.37 * 1e-7 ,2.5 * 1e-7, 2.62 * 1e-7, 2.66 * 1e-7];
vg08_qc_10nm=[0,1.51 * 1e-7, 1.8 * 1e-7, 1.9 * 1e-7, 2.0 * 1e-7, 2.05 * 1e-7];

vg10_noqc_6nm=[0, 2.94 * 1e-7, 3.4 * 1e-7,3.59 * 1e-7, 3.7 * 1e-7, 3.7 * 1e-7];
vg10_qc_6nm=[0, 1.85 * 1e-7, 2.2 * 1e-7,0,2.5 * 1e-7, 2.5 * 1e-7];

%vg1_noqc_6nm=[2.94 * 1e-7, 3.4 * 1e-7,3.59 * 1e-7, 3.7 * 1e-7, 3.7 * 1e-7];

%vg0.8_noqc_6nm=[1.74 * 1e-7, 1.92 * 1e-7, 0, 2.07 * 1e-7, 2.66 * 1e-7];

vg10_noqc_10nm = 100 * vg10_noqc_10nm;
vg10_qc_10nm = 100 * vg10_qc_10nm;
vg10_noqc_6nm = 100 * vg10_noqc_6nm;
vg10_qc_6nm = 100 * vg10_qc_6nm;

hold on
h(1)=plot(vd, vg10_noqc_10nm,'k^-','LineWidth', 6, 'MarkerSize', 12);
h(2)=plot(vd, vg10_qc_10nm,'k^--','LineWidth',6, 'MarkerSize', 12);
h(3)=plot(vd, vg10_noqc_6nm,'ko-','LineWidth', 6, 'MarkerSize', 12);
h(4)=plot(vd, vg10_qc_6nm,'ko--','LineWidth',6, 'MarkerSize', 12);

set(gca,'FontSize',25);
set(gca,'LineWidth',3);

%xlabel('V_d / V','FontSize', 30);
ylabel('Drain Current (A)','FontSize', 30);

xlim([0 1.1]);
ylim([0 6e-5]);
set(gca,'XTick', [0, 0.2, 0.4, 0.6, 0.8, 1.0]);

text(0.35, -7.0e-6, 'Drain Voltage V_d (V)' , 'FontSize', 30);
%text(0.4,1e-7,'V_g=0.6 V, noqc', 'FontSize', 30);
%legend('V_g=0.6 V, noqc','V_g=0.8 V, noqc', 'V_g=1 V, noqc','V_g=0.6 V, qc', 'V_g=0.8 V, qc', 'V_g=1 V, qc');
[legh, objh, outh, outm] = legend('W=10nm, noqc','W=10nm, qc','W=6nm, noqc', 'W=6nm, qc');
%set(legh,'Box','off');
%legh2=copyobj(legh,gcf);
%[legh2,objh2]=legend(h(4:6),'V_g=0.6 V, qc', 'V_g=0.8 V, qc', 'V_g=1 V, qc',2);
%set(legh2,'Box','off');
