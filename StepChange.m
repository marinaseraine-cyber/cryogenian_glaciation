subplot(2,1,1)
a = [median(detd18o) 0.1 1 mad(diff(detd18o),1)*1];

[M,P,~] = CPRBayes(detd18o,'normal','alpha',a,'conf',0.95,'depth',1,'paramtype','bayesian');

d_fit = zeros(length(detd18o),3);
for i = 1:length(M)-1
    d_fit(M(i)+1:M(i+1),1) = repmat(P(i,1),M(i+1)-M(i),1);
    d_fit(M(i)+1:M(i+1),2) = d_fit(M(i)+1:M(i+1),1) + P(i,2);
    d_fit(M(i)+1:M(i+1),3) = d_fit(M(i)+1:M(i+1),1) - P(i,2);
end

scatter(detages,detd18o,5,'o','k','filled'); hold on;
%errorbar(Age,detd18o,RbF_sigma, 'LineStyle','none','Color','k'); hold on
plot(detages,d_fit(:,1),'b'); hold on;
plot(detages,d_fit(:,2),'r'); hold on;
plot(detages,d_fit(:,3),'r'); hold off;
axis([500 900 -60 40])
set(gca, 'XDir','reverse')
xlabel('Age (Ma)') 
ylabel(['\delta^{18}O (',char(8240),')'])
% x = [0.73 0.78];
% y = [0.45 0.41];
% annotation('textarrow',x,y,'String','2.1 Ga')
legend('data','CPR median','CPR 95% c.i.','Location', 'Northwest' );
title('X \delta^{18}O (Database)')

% subplot(2,1,2)
% a = [median(detd18o2) 0.1 1 mad(diff(detd18o2),1)*1];
% 
% [M,P,~] = CPRBayes(detd18o2,'normal','alpha',a,'conf',0.95,'depth',5,'paramtype','bayesian');
% 
% d_fit = zeros(length(detd18o2),3);
% for i = 1:length(M)-1
%     d_fit(M(i)+1:M(i+1),1) = repmat(P(i,1),M(i+1)-M(i),1);
%     d_fit(M(i)+1:M(i+1),2) = d_fit(M(i)+1:M(i+1),1) + P(i,2);
%     d_fit(M(i)+1:M(i+1),3) = d_fit(M(i)+1:M(i+1),1) - P(i,2);
% end
% 
% scatter(detages2,detd18o2,5,'o','k','filled'); hold on;
% %errorbar(Age,detd18o,RbF_sigma, 'LineStyle','none','Color','k'); hold on
% plot(detages2,d_fit(:,1),'b'); hold on;
% plot(detages2,d_fit(:,2),'r'); hold on;
% plot(detages2,d_fit(:,3),'r'); hold off;
% axis([500 900 0 15])
% set(gca, 'XDir','reverse')
% xlabel('Age (Ma)') 
% ylabel(['\delta^{18}O (',char(8240),')'])
% % x = [0.73 0.78];
% % y = [0.45 0.41];
% % annotation('textarrow',x,y,'String','X Ga')
% legend('data','CPR median','CPR 95% c.i.','Location', 'Northwest' );
% title('X \delta^{18}O (Database)')