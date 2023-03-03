load data_customer.txt;
data=data_customer(:,1:2); % 提取出出发的数据
num=length(data); 
aini=0.05;
sigmaini=0.005;
LO=sqrt(10)/2;
T=400;  % 迭代次数 可改小 大概300次拟合
t=1; 
i=1;
for j1=1:5
    for j2=1:2
        w1(j1,j2)=113.4+(114.6-113.4)*rand;
        w2(j1,j2)=22.4+(23.4-22.4)*rand;
        i=i+1;
    end
end
x1 = data(:,1);
x2 = data(:,2);
while (t<=T)
    a=aini*exp(-t/T);
    sigma=sigmaini*exp(-t/T);
    LN=round(LO*(1-t/T));
    for i=1:num
        e_norm=(x1(i)-w1).^2+(x2(i)-w2).^2;
        minj1=1;minj2=1;
        min_norm=e_norm(minj1,minj2);
        for j1=1:5
            for j2=1:2
                if e_norm(j1,j2)<min_norm
                    min_norm=e_norm(j1,j2);
                    minj1=j1;
                    minj2=j2;
                end
            end
        end
        j1_c= minj1;
        j2_c= minj2;
        bias=[1 2;3 4; 5 6; 7 8;9 10];
        idx(i)=bias(j1_c,j2_c);
        w1(j1_c,j2_c)=w1(j1_c,j2_c) + a * (x1(i) - w1(j1_c,j2_c));
        w2(j1_c,j2_c)=w2(j1_c,j2_c) + a * (x2(i) - w2(j1_c,j2_c));
        for neighbour_radius=1:1:LN
            jj1=j1_c - neighbour_radius;
            jj2=j2_c;
            if ((jj1>=1)&&(jj1<=5))
                e_factor = exp(-((jj1-j1_c).^2+(jj2-j2_c).^2)/2*sigma^2);
                w1(jj1,jj2)=w1(jj1,jj2) + a * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + a * e_factor * (x2(i)-w2(jj1,jj2));
            end
            jj1=j1_c + neighbour_radius;
            jj2=j2_c;
            if ((jj1>=1)&&(jj1<=5))
                e_factor = exp(-((jj1-j1_c).^2+(jj2-j2_c).^2)/2*sigma^2);
                w1(jj1,jj2)=w1(jj1,jj2) + a * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + a * e_factor * (x2(i)-w2(jj1,jj2));
            end
            jj1=j1_c;
            jj2=j2_c - neighbour_radius;
            if ((jj2>=1)&&(jj2<=2))
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma^2);
                w1(jj1,jj2)=w1(jj1,jj2) + a * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + a * e_factor * (x2(i)-w2(jj1,jj2));
            end
            jj1=j1_c;
            jj2=j2_c + neighbour_radius;
            if ((jj2>=1)&&(jj2<=2))
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma^2);
                w1(jj1,jj2)=w1(jj1,jj2) + a * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + a * e_factor * (x2(i)-w2(jj1,jj2));
            end
        end
    end
    t=t+1;
    figure(1)
    plot(x1,x2,'.','MarkerSize',10)
    hold on
    plot(w1,w2,'.','MarkerSize',30,'Color',[0.4 0.4 0.6])
    plot(w1,w2,'b','linewidth',1)
    plot(w1',w2','b','linewidth',1)
    hold off
    title(['t=' num2str(t),' times']);
    drawnow
end
hold on
c =[1 0.5256 0;1 0 1; 0 1 1;1 0 0; 0 1 0;0 0 1;0 0 0;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
gscatter(data(:,1),data(:,2),idx,c);
