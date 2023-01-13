 % Roofline Figure
[figSettings, color, rgbcmy, mList, lineList, dashList, subplot_pos] = setFigureDefaultsNew();

% for a serial benchmark on the hardware
peakflops = 5546.22;
peakmem = 2773.11;

x = 0:0.01:10;

roof = x*peakmem;
roof(roof>peakflops) = peakflops;

darkgreen = [0.1797 0.5430 0.3398];

flops_paper = 441.7;
intensity_paper = 0.163;

flops_15M = 456.8;
intensity_15M = 0.173;

flops_8M = 426.5;
intensity_8M = 0.159;

flops_6M = 412.2;
intensity_6M = 0.155;

flops_4M = 373.6;
intensity_4M = 0.145;

flops_2M = 308.7;
intensity_2M = 0.117;

% double particle density with same sized domain
flops_bigdensity = 433.7;
intensity_bigdensity = 0.242;

intensity = [intensity_bigdensity intensity_2M intensity_4M intensity_6M intensity_8M intensity_15M];
flops = [flops_bigdensity flops_2M flops_4M flops_6M flops_8M flops_15M];

h = plot(x, roof, 'k-');
hold on
plot(intensity_paper, flops_paper, '.', 'MarkerSize', 30)
scatter(intensity, flops, 75, darkgreen,'+')
set(h, 'MarkerFaceColor', get(h,'Color'));
xline(peakflops/peakmem,'--','LineWidth',3)
grid on
xlabel('Arithmetic Intensity (FLOPs/Byte)')
ylabel('Performance (MFLOPs/s)')
set(gca,'xscale','log')
set(gca,'yscale','log')
rectangle('Position', [0.1 300 0.15 200])
ylim([10^2 10^4])
legend('Roofline Model','Parameter Set for Speedup Results','Varying Particle Number','Location','northwest')
ax1 = gca;

ax2 = axes('Position',[.58 .2 .3 .27]);
box on;
plot(x, roof, 'k-');
hold on
plot(intensity_paper, flops_paper, '.', 'MarkerSize', 45)
scatter(intensity, flops, 125, darkgreen,'+')
grid on;
set(ax2,'xscale','log')
set(ax2,'yscale','log')
xlim([0.1 0.25])
ylim([300 500])
set(ax2,'xticklabel',[],'yticklabel',[])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% saveas(gcf,'RooflineInset','epsc')