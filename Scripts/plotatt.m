figure
plot(record_att);
xlabel('time');
ylabel('attention level');
yy=mean(record_att)
l=length(record_att);
for i=1:l
   y(i)=yy; 
end
hold on
plot(y,'r')
ylim([0, 10]);
legend('real-time attention level','average attention level')
saveas(gcf,'att.jpg')
figure
plot(record_med);
xlabel('time');
ylabel('meditation level');
yym=mean(record_med)
l=length(record_med);
for i=1:l
   ym(i)=yym; 
end
hold on
plot(ym,'r')
ylim([0, 10]);
legend('real-time meditation level','average meditation level')
saveas(gcf,'med.jpg')
