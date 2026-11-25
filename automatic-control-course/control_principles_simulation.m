function control_principles_simulation(showFigures)
% Reproduce Figures 2-7 and Tables 2-3 of the revised paper.
% In MATLAB, change to this file's folder and run:
%   control_principles_simulation
% For a numerical run without figures:
%   control_principles_simulation(false)
% Uses base MATLAB only. Time is in s and angular frequency in rad/s.

if nargin == 0
    showFigures = true;
end

%% 4.1 Modeling: unity negative feedback, K = 20
openDen = [1 7 14 8];
seriesDen = conv(conv([1 1], [1 2]), [1 4]);
assert(isequal(seriesDen, openDen));
K = 20;
closedDen = openDen + [0 0 0 K];
seriesClosedDen = seriesDen + [0 0 0 K];
assert(isequal(seriesClosedDen, closedDen));
fprintf('L(s) numerator = %g, denominator = [%s]\n', K, num2str(openDen));
fprintf('T(s) numerator = %g, denominator = [%s]\n', K, num2str(closedDen));
disp('Closed-loop poles at K = 20:'); disp(roots(closedDen));
exerciseDen = conv(conv([1 1], [1 2]), [1 3]) + [0 0 0 K];
fprintf('Exercise: move -4 to -3; closed-loop denominator = [%s]\n', num2str(exerciseDen));

if showFigures
    figure('Name', 'Figure 2 - Unity negative feedback', 'NumberTitle', 'off');
    axes('Position',[.05 .12 .9 .75]); hold on; axis([0 10 0 4]); axis off;
    rectangle('Position',[1.35 2.15 .7 .7], 'Curvature',[1 1]);
    text(1.7,2.5,'+','HorizontalAlignment','center');
    text(1.75,2.02,'-','HorizontalAlignment','center');
    rectangle('Position',[2.65 2.05 1.15 .9]);
    text(3.23,2.5,'K','HorizontalAlignment','center');
    rectangle('Position',[4.45 2.05 3.5 .9]);
    text(6.2,2.5,'1 / [(s+1)(s+2)(s+4)]','HorizontalAlignment','center');
    quiver([.45 2.05 3.8 7.95],[2.5 2.5 2.5 2.5],[.9 .6 .65 1.55],zeros(1,4),0,'k');
    plot([8.6 8.6 1.7],[2.5 .85 .85],'k');
    quiver(1.7,.85,0,1.3,0,'k');
    plot(8.6,2.5,'k.','MarkerSize',14);
    text(.4,2.87,'R(s)'); text(8.9,2.87,'Y(s)');
    text(4.9,.55,'H(s) = 1','HorizontalAlignment','center');
end

%% 4.2 Second-order step response and performance indices
wn = 2;
zeta = [0.2 0.4 0.7 1.0 1.5];
t = (0:0.001:40).';
y = zeros(numel(t), numel(zeta));
overshoot = zeros(size(zeta));
peakTime = NaN(size(zeta));  % No finite peak for monotone responses
settlingTime = zeros(size(zeta));
for i = 1:numel(zeta)
    z = zeta(i);
    if z < 1
        wd = wn * sqrt(1-z^2);
        y(:,i) = 1 - exp(-z*wn*t) .* (cos(wd*t) + z*wn/wd*sin(wd*t));
        peakTime(i) = pi/wd;
        overshoot(i) = 100*exp(-z*wn*peakTime(i));
    elseif z == 1
        y(:,i) = 1 - exp(-wn*t) .* (1+wn*t);
    else
        p1 = -wn/(z+sqrt(z^2-1));
        p2 = -wn*(z+sqrt(z^2-1));
        y(:,i) = 1 + (p2*exp(p1*t)-p1*exp(p2*t))/(p1-p2);
    end
    lastOutside = find(abs(y(:,i)-1)>0.02,1,'last');
    assert(lastOutside < numel(t), 'Increase the simulation horizon.');
    settlingTime(i) = t(lastOutside+1);
end
fprintf('\nTable 2: damping, overshoot (%%), peak time (s), 2%% settling time (s)\n');
fprintf('%5.1f %10.2f %12.3f %12.3f\n', [zeta; overshoot; peakTime; settlingTime]);
if showFigures
    figure('Name','Figure 3 - Second-order step response','NumberTitle','off');
    plot(t,y,'LineWidth',1.3); grid on; xlim([0 12]);
    xlabel('Time (s)'); ylabel('Unit-step response');
    labels = arrayfun(@(z) sprintf('\\zeta = %.1f',z), zeta,'UniformOutput',false);
    legend(labels,'Location','northeast');
end

%% 4.3 Root locus and third-order closed-loop step response
gains = [20 60 90 110];
colours = lines(numel(gains));
criticalGain = 7*14-8;
criticalFrequency = sqrt(14);
t3 = (0:0.01:18).';
y3 = zeros(numel(t3),numel(gains));
for i = 1:numel(gains)
    K = gains(i);
    den = openDen + [0 0 0 K];
    p = roots(den);
    % Step output transform K / [s*D(s)]. Poles are simple for these gains.
    % Residue at s=0 gives the constant term, even in unstable cases.
    coefficients = K ./ (p .* polyval(polyder(den),p));
    y3(:,i) = K/den(end) + real(exp(t3*p.')*coefficients);
    fprintf('\nK = %g; closed-loop poles:\n',K); disp(p);
end
if showFigures
    scan = unique([0 linspace(0,160,1601) logspace(log10(160),3.5,700) criticalGain]);
    locus = zeros(3,numel(scan));
    for i = 1:numel(scan)
        locus(:,i) = roots(openDen + [0 0 0 scan(i)]);
    end
    figure('Name','Figure 4 - Root locus','NumberTitle','off');
    plot(real(locus(:)),imag(locus(:)),'.','MarkerSize',5); hold on; grid on;
    h0 = plot([-1 -2 -4],[0 0 0],'kx','MarkerSize',9,'LineWidth',1.5);
    hc = plot([0 0],[-1 1]*criticalFrequency,'ro','MarkerSize',8,'LineWidth',1.5);
    limY = [-8 8]; plot([0 0],limY,'k:'); ylim(limY); xlim([-11 4]);
    xlabel('Real part'); ylabel('Imaginary part');
    legend([h0 hc],{'Open-loop poles','K = 90 crossing'},'Location','northwest');

    figure('Name','Figure 5 - Third-order step response','NumberTitle','off');
    set(gca,'ColorOrder',colours,'NextPlot','replacechildren');
    plot(t3,y3,'LineWidth',1.3); grid on;
    xlabel('Time (s)'); ylabel('Unit-step response');
    gainLabels = arrayfun(@(k) sprintf('K = %g',k),gains,'UniformOutput',false);
    legend(gainLabels,'Location','northwest');
end

%% 4.4 Frequency responses and stability margins
w = unique([logspace(-2,2,1200), criticalFrequency]);
frequencyDen = polyval(openDen,1i*w);
Krows = gains(:);
response = Krows*(1./frequencyDen);
phase = -(atan(w)+atan(w/2)+atan(w/4))*180/pi;
gainDB = 20*log10(criticalGain./gains);
crossing = zeros(size(gains));
phaseMargin = zeros(size(gains));
for i = 1:numel(gains)
    K = gains(i);
    crossing(i) = fzero(@(x) abs(polyval(openDen,1i*x))-K,[0 100]);
    wc = crossing(i);
    phaseMargin(i) = 180-(atan(wc)+atan(wc/2)+atan(wc/4))*180/pi;
end
phaseMargin(abs(phaseMargin)<1e-9) = 0;
fprintf('\nTable 3: K, gain crossover (rad/s), phase margin (deg), gain margin (dB)\n');
fprintf('%5.0f %12.3f %12.2f %12.2f\n', [gains; crossing; phaseMargin; gainDB]);
qualifies = phaseMargin>=45 & gainDB>=6;
fprintf('Gain(s) satisfying both margin requirements: %s\n',num2str(gains(qualifies)));
fprintf('Stable step final value at qualifying gain(s): %s\n', ...
    num2str(gains(qualifies)./(8+gains(qualifies))));

if showFigures
    figure('Name','Figure 6 - Bode plots','NumberTitle','off');
    subplot(2,1,1); hold on; grid on;
    set(gca,'XScale','log');
    for i = 1:numel(gains)
        semilogx(w,20*log10(abs(response(i,:))),'Color',colours(i,:),'LineWidth',1.3);
    end
    legend(gainLabels,'Location','southwest'); ylabel('Magnitude (dB)');
    plot([w(1) w(end)],[0 0],'k:','HandleVisibility','off');
    xlim([w(1) w(end)]);
    subplot(2,1,2);
    semilogx(w,phase,'k','LineWidth',1.3); grid on; hold on;
    plot([w(1) w(end)],[-180 -180],'k:');
    xlim([w(1) w(end)]);
    xlabel('Angular frequency (rad/s)'); ylabel('Phase (deg)');

    % For this strictly proper loop, the infinite-radius arc maps to zero.
    % Include DC and the zero limit when tracing the full Nyquist contour.
    wNyq = unique([0 logspace(-4,4,3200) criticalFrequency]);
    figure('Name','Figure 7 - Nyquist detail and full curve','NumberTitle','off');
    for panel = 1:2
        subplot(1,2,panel); hold on; grid on;
        for i = 1:numel(gains)
            positive = gains(i)./polyval(openDen,1i*wNyq);
            contour = [0 conj(fliplr(positive(2:end))) positive 0];
            plot(real(contour),imag(contour),'Color',colours(i,:),'LineWidth',1.3);
        end
        plot(-1,0,'kx','MarkerSize',9,'LineWidth',1.5,'HandleVisibility','off');
        xlabel('Real part'); ylabel('Imaginary part');
        if panel == 1
            xlim([-1.9 -0.45]); ylim([-0.75 0.75]); title('Critical-point detail');
        else
            axis equal; title('Full Nyquist curve');
            legend(gainLabels,'Location','northeast');
        end
    end
end
end
