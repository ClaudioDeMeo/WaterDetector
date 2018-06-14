% insert all shapefiles in a directory together this script
files=dir('*.shp'); 
total=zeros(length(files),1);
date=string(total);

for i = 1:length(files)
    total(i)=length(shaperead(files(i).name));
    date(i)=extractBefore(files(i).name,9);
end 

bar(categorical(date),total);
hold on;
med=ones(length(files),1) * mean(total);

p=plot(categorical(date),med);
p.LineWidth=2;
title('Tanks number variation');
xlabel('date');
ylabel('number of tanks');
legend('number of tanks',strcat('average number of tanks: ',string(round(mean(total)))));
hold off;