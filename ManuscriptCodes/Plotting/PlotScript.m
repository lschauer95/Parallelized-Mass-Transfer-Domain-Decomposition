[figSettings, color, rgbcmy, mList, lineList, dashList, subplot_pos] = setFigureDefaultsNew();

load('PlotData.mat')

    %=============================================================%
    %===================== FIGURE SCRIPTS ========================%

                        %Figure 3 (error figure)

%         h = plot(smallerrcores, smallerrors, '-*');
%         ylim([1.5*10^(-3) 2.5*10^(-3)])
%         grid on
%         xlabel('Number of Cores')
%         ylabel('RMSE')
%         set(gca,'yscale','log')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['Error2DHeavisideCheck10M'],'epsc')


%                           %Figure 6(a) 
%         h = plot(Ns2D10M, MTwork2D10M, '-d', Ns2D10M, alphafunc2D10M, '-');
%         set(gca,'FontSize',34)
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         xlabel('$N_S$')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^3, 10^6], 'ylim', [10^(-1), 10^3]);
%         legend('Observed MT Time 2-d','$\alpha_2N_s\log(N_s)$','Location','northwest','Interpreter','latex')
%         dim = [.65 .25 .3 .3];
%         str = '\alpha_2 = 1.054e-06';
%         annotation('textbox',dim,'String',str,'FitBoxToText','on');
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['NsvsMT2D10M'],'epsc')


%                    Figure 6(b)         

%         h = plot(Ns3D5M, MTwork3D5M, '-d', Ns3D5M, alphafunc3D5M, '-');
%         set(gca,'FontSize',34)
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         xlabel('$N_S$')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^3, 10^6], 'ylim', [10^(-1), 10^3]);
%         dim = [.65 .25 .3 .3];
%         str = '\alpha_3 = 1.707e^{-4}';
%         annotation('textbox',dim,'String',str,'FitBoxToText','on');
%         legend('Observed MT Time 3-d','$\alpha_3N_s\log(N_s)$','Location','northwest','Interpreter','latex')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['NsvsMT3D5M'],'epsc')



%                   %FIGURE 7(a)
         
%         h = plot(Ns2D10M, MTwork2D10M, '-d',...
%                 Ns2D15M, MTwork2D15M, '-o',...
%                 Ns2D20M, MTwork2D20M, '-s');
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         xlabel('$N_S$')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^3, 10^6], 'ylim', [10^(-1), 10^3]);
%         legend('$N=10M$ 2-d','$N=15M$ 2-d','$N=20M$ 2-d','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['NsvsMT2Dfull'],'epsc')


%                   %FIGURE 7(b)

%         h = plot(cores, MTwork2D10M, '-d', cores, 10^4*1./cores, '-');
%         grid on
%         xlabel('Number of Cores')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^1, 10^4], 'ylim', [1, 10^3]);
%         legend('$N=10M$ 2-d','Linear Speedup','Location','northeast')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['strongscaling2d'],'epsc')
    

%                     %Figure 8(a)
% 
%         h = plot(Ns3D5M, MTwork3D5M, '-d', ...
%             Ns3D10M, MTwork3D10M, '-o', Ns3D15M, MTwork3D15M, '-s');
% 
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));       
%         grid on
%         xlabel('$N_S$')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^3, 10^6], 'ylim', [10^(-1), 10^3]);
%         legend('$N=5M$ 3-d','$N=10M$ 3-d','$N=15M$ 3-d','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['NsvsMT3Dfull'],'epsc')


%                   %Figure 8(b)
%
%         h = plot(cores3D, MTwork3D5M, '-d', cores3D, 10^4*1./cores3D, '-');
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         xlabel('Number of Cores')
%         ylabel('MT Time (s)')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^1, 10^4], 'ylim', [10^(0), 10^3]);
%         legend('$N=5M$ 3-d','Linear Speedup', 'Location','northeast')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['strongscaling3d'],'epsc')


%                   %Figure 9

%     h = plot(lowcores(2:4), lowcoresMPItime(2:4), '-d');
%     grid on
%     xlabel('Number of Cores')
%     ylabel('MPI Time (s)')
%     xticks(lowcores(2:4))
%     set(gca,'xlim', [0, 20]);
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %     saveas(gcf, ['singlenodeMPIscaling'],'epsc')


                        %Figure 10(a)
                        
%         h = plot(Ns2D10M, MPIwork2D10M, '-x', Ns2D15M, MPIwork2D15M, '-*', ...
%             Ns2D20M, MPIwork2D20M, '-+', Ns3D5M, MPIwork3D5M, '-d', ...
%             Ns3D10M, MPIwork3D10M, '-o', Ns3D15M, MPIwork3D15M, '-s');
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));       
%         grid on
%         xlabel('$N_S$')
%         ylabel('MPI Time (s)')
%         legend('$N=10M$ 2-d','$N=15M$ 2-d','$N=20M$ 2-d','$N=5M$ 3-d','$N=10M$ 3-d',...
%             '$N=15M$ 3-d','Location','northwest')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^3, 10^6], 'ylim', [10^(-1), 10^3]);
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['NsvsSwapallrunsloglog'],'epsc')


%                    %Figure 10(b)

%         h = plot(cores, MPIwork2D10M, '-d', cores3D, MPIwork3D5M, '-o');
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         xlabel('Number of Cores')
%         ylabel('MPI Time (s)')
%         legend('$N=10M$ 2-d','$N=5M$ 3-d','Location','northwest')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         set(gca,'xlim', [10^1, 10^4], 'ylim', [10^(-1), 10^2]);
%         text(10^4.1,10^2.7,'$\mathbf{\leftarrow}$ Increasing cores')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['MPIstrongscaling'],'epsc')


                    % Figure 11(a)

%         P=1:10:5001;
%         L1 = 1000000;
%         L2 = sqrt(1000000);
%         L3 = (1000000)^(1/3);
%         S1D=(1./P + 2*xi/L1).^(-1);
%         S2D=(P.^(-0.5) + 2*xi/L2).^(-2);
%         S3D=(P.^(-0.3333) + 2*xi/L3).^(-3);
% 
%         h = plot(P, S1D, '-', P, S2D, '-', P, S3D, '-', P, P, 'k-');
%         hold on 
%         xline(1773,'--','Linewidth',3.5,'Color', [0.85 0.325 0.098])
%         xline(321,'--','Linewidth',3.5,'Color', [0.929 0.694 0.125])
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         xlim([0 5000])
%         legend('Projected 1-d Speedup','Projected 2-d Speedup','Projected 3-d Speedup','Perfect 1:1 Speedup',...
%             '2-d Projection for $\mathcal{E}=0.75$','3-d Projection for $\mathcal{E}=0.5$','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['ProjectedSpeeduphyper'],'epsc')


%               % Figure 11(b)

%         P = 1:10:5001;
%         L1 = 1000;
%         L2 = 1000;
%         L3 = 1000;
%         S1D=(1./P + 2*xi/L1).^(-1);
%         S2D=(P.^(-0.5) + 2*xi/L2).^(-2);
%         S3D=(P.^(-0.3333) + 2*xi/L3).^(-3);       
% 
%         h = plot(P, S1D, '-', P, S2D, '-', P, S3D, '-', P, P, 'k-');
%         hold on 
%         xline(1773,'--','Linewidth',3.5,'Color', [0.85 0.325 0.098])
%         xline(3156,'--','Linewidth',3.5,'Color', [0.929 0.694 0.125])
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         xlim([0 5000])
%         legend('Projected 1-d Speedup','Projected 2-d Speedup','Projected 3-d Speedup','Perfect 1:1 Speedup',...
%             '2-d Projection for $\mathcal{E}=0.75$','3-d Projection for $\mathcal{E}=0.5$','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['ProjectedSpeeduphyper'],'epsc')


                        %Figure 12
                        
%         %non-normalized speedup using serial time t=15819
%         P=1:10:2701;
%         S2D=(P.^(-0.5) + 2*xi/L2).^(-2);
% 
%         h = plot(cores, speedup2D10Mserial, '-d', P, S2D, '-', cores, cores, 'k-');
%         hold on
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xline(1773,'--','Linewidth',4)
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         legend('Checkerboard Speedup','Projected 2-d Speedup','1:1','Efficiency Condition for $\mathcal{E}=0.75$','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['speedup2D10Mserial'],'epsc')


                        %Figure 13
%         d = 3;
%         eff = 0.5;
%         eta = 0.5*((1/eff)^(1/d) - 1);
%         wdom = xi/eta;
%         P=1:10:1731;
%         S3D=(P.^(-0.3333) + 2*xi/L3).^(-3);
%         projcores = (100/wdom)^d;
%                          
%         h = plot(fullcores3D, fullspeedup3D5M, '-d', P, S3D, '-', fullcores3D, fullcores3D, 'k-');
%         hold on 
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xline(projcores,'--','Linewidth',4)
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         legend('Speedup 3-d 5M','Projected 3-d Speedup','1:1','Efficiency Condition for $\mathcal{E}=0.5$','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['full3D5Mspeedup'],'epsc')



                        %Figure 13

%         L1 = 1000;
%         P=1:10:501;
%         S1D=(1./P + 2*xi/L1).^(-1);
% 
%         h = plot(slicescores(1:9), slicesspeedup(1:9), '-d', P, S1D, '-', slicescores(1:9), slicescores(1:9), 'k-');
%         set(h, {'MarkerFaceColor'}, get(h,'Color'));
%         grid on
%         hold on
%         xline(89, '--', 'Linewidth', 4)
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         xlim([0 500])
%         legend('Slices Speedup','Projected Slices Speedup','1:1','Efficiency Condition for $\mathcal{E}=0.75$','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['slicesspeedup10M'],'epsc')



                        %Figure 15

%         P = 1:10:1601;
%         L1 = 1000;
%         L2 = 1000;
%         S1D=(1./P + 2*xi/L1).^(-1);
%         S2D=(P.^(-0.5) + 2*xi/L2).^(-2);
%         
%         h = plot(badcores, badcheckspeedup, '-d', P, S1D, '-', P, S2D, '-', badcores, badcores, 'k-');
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         xticks(badcores)
%         xlim([0 1600])
%         ylim([0 1600])
%         legend('Checkerboard Speedup','1:1','Predicted 1-d Speedup','Predicted 2-d Speedup','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['Speedup2DHeavisideCheck10Mfullbadexample'],'epsc')


                        %Figure 16
                        
%         h = plot(cores, speedup2D10M, '-d', cores, speedup2D15M, '-o',...
%              cores, speedup2D20M, '-s', cores, speedupcores2D, 'k-');
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         legend('$N=10M$ 2-d','$N=15M$ 2-d','$N=20M$ 2-d','1:1','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['speedup2Drelativeall'],'epsc')



                    %Figure 17
                    
%         h = plot(cores3D, speedup3D5M, '-d', cores3D, speedup3D10M, '-o',...
%              cores3D15M, speedup3D15M, '-s', cores3D, speedupcores3D, 'k-');
%         set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%         grid on
%         xlabel('Number of Cores')
%         ylabel('Speedup')
%         legend('$N=5M$ 3-d','$N=10M$ 3-d','$N=15M$ 3-d','1:1','Location','northwest')
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% %         saveas(gcf, ['speedup3Drelativeall'],'epsc')