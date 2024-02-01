%People generation with borders as any equation. Speed fluctuates depending
%on the Area (higher speed for lower area) and distance between nodes
%fluctuates dependding on Area (higher distance for higher area). The high
%speed is associated with the nodes running from the walls and the low
%speed with the nodes in region II`

clear all;
close all

%%

Xmin=-2;           %range of image
Ymin=-4;
Xmax=20;
Ymax=3;

Nk=40;      % no. of nodes
dt=0.05;     %time increment / time step
time=200; %Number of iterations

R=15;   % neighborhood radius
rmax=2;  % neighbors' distance from each other
rmin=1;


lamda=0.5;  %people speed and movement parameter/ variable eq. final velocity of people
alpha=1;  % people speed and movement parameter
gamma=1;  % people speed and movement parameter /variable eq. final velocity of people
mu=0.5;  %size of step nodes can move in one interval/variable eq. 23

COFv=0; %coefficient of the speed
COFvmin=0.6; %speed of nodes
COFvmax=3.6; %speed of nodes
Vo=2.1; %original speed
Vnew=0;


%distance between k and target, used to make some nodes go faster than
%others (tail effect)
COFdistmax=0.1; %min
COFdistmin=4; %max

%Max and Min chord
Lcmin=1;
Lcmax=10;

%% BORDERS 
syms Z1 Z2
f11= 2*cos(((2*pi)/30)*(Z1-1.5))+1;
f22= 4*sin(((2*pi)/40)*(Z2+4))-3.5;

f1= @(Z)2*cos(((2*pi)/30)*(Z-1.5))+1;
f2= @(Z)4*sin(((2*pi)/40)*(Z+4))-3.5;

%% To calculate the chords
Step=0.1;
xCoordinate=Xmin-5;
for i=1:((Xmax-(Xmin-5))/Step) 
        xCoordinate=xCoordinate+Step;
        PointA(1,i)=real(xCoordinate);
        PointA(2,i)=real(subs(f11, Z1, xCoordinate));
        PointB(1,i)=real(xCoordinate);
        PointB(2,i)=real(subs(f22, Z2, xCoordinate));
        MidPoints(:,i) = real((round(PointA(:,i),2) + round(PointB(:,i),2)).'/2);

end

%%
rnoise=cell(1,time);
COFr=cell(1,time);
va=cell(1,time+1);
vb=cell(1,time+1);
delta=cell(1,time);
v=cell(1,time+1);
v{1,1}=zeros(2,Nk);
phi=cell(1,time+1);
x{1,1}=randn(2,Nk);


%%
for i=1:Nk
    xCoord=x{1,1}(1,i);
    yCoord=x{1,1}(2,i);
    if (yCoord >= f1(xCoord))
        x{1,1}(2,i)=f1(xCoord)-1;
    end
    if (yCoord <= f2(xCoord))
        x{1,1}(2,i)=f2(xCoord)+1;
    end
end

%%


vg=cell(1,time+1);
vg{1,1}=zeros(2,Nk);
Lck=zeros(time, Nk);



wt1=[6;1];
wt2=[13;-1];
wt3=[18; -3];
wt=wt1;

%% program for nodes (Main program)
for i=1:time
   
    vg{1,i+1}=zeros(2,Nk);
    delta{1,i}=zeros(2,Nk);
    rnoise{1,i}=zeros(2,Nk);
    COFr{1,i}=zeros(2,Nk);
    
     %analise distance between k and target. K1 should go faster than Kn
     %for each time i, the distance for all k     
    for k=1:Nk    
     dist1(i,k)=norm(wt-x{1,i}(:,k));
    end
 
           
    for k=1:Nk
%calculate chord through node xk at time i
        for j=1:length(MidPoints)
            a=MidPoints(:,j);
            if (a(1)==round(x{1,i}(1,k),1))
                PointC=MidPoints(:,j);
                PointD=MidPoints(:,j-1);
                coefficients = polyfit(PointC, PointD, 1); %eq between Point C and D
                fPerpendicularMidPointAB1=(-1/coefficients(1))*(Z1-x{1,i}(1,k))+x{1,i}(2,k);
                fPerpendicularMidPointAB2=(-1/coefficients(1))*(Z2-x{1,i}(1,k))+x{1,i}(2,k);
                x1=double(vpasolve(fPerpendicularMidPointAB1==f11,Z1));
                x2=double(vpasolve(fPerpendicularMidPointAB2==f22,Z2));
                y1=double(subs(f11,Z1,x1));
                y2=double(subs(f22,Z2,x2));
                Point1=[x1;y1];
                Point2=[x2;y2];
                Lck(i,k)= double(norm (Point1-Point2)); %chrod of k at time i
                if (Lck(i,k)<Lcmin)
                    Lck(i,k)=Lcmin;
                elseif (Lck(i,k)>Lcmax)
                    Lck(i,k)=Lcmax;
                end
            end
        end
        
        Vnew(i,k)=((Lcmax - Lck(i,k))/(Lcmax-Lcmin))*(COFvmax-COFvmin)+COFvmin;

        if (k==5) %half the speed for the disabled node
            Vnew(i,k)=Vnew(i,k)/2;
        end

         Z=x{1,i}(1,k);
         a1=[Z;(f1(Z)-0.5)];
         a2=[Z;(f2(Z)+0.5)];
 
        if (x{1,i}(2,k)<(f1(Z)) && x{1,i}(2,k)>(f2(Z))) %region II %we evaluate the values of y of the node and the border at the same x
            va{1,i+1}(:,k)=Vnew(i,k)*((wt-x{1,i}(:,k))./norm(wt-x{1,i}(:,k)));
            
        else 
               SS1 = norm(a1-x{1,i}(:,k));
               if SS1>3
                SS1=3; 
               end
               if SS1<1
                SS1=1; 
               end
               SS2 = norm(a2-x{1,i}(:,k));
               if SS2>3
                SS2=3; 
               end
               if SS2<1
                SS2=1; 
               end

            if (x{1,i}(2,k)>=(f1(Z))) %region I

                va{1,i+1}(:,k)=-Vnew(i,k)*((R/SS1)-1)*(x{1,i}(:,k)-a1);
                
            elseif (x{1,i}(2,k)<=(f2(Z))) %region III

                va{1,i+1}(:,k)=-Vnew(i,k)*((R/SS2)-1)*(x{1,i}(:,k)-a2);

            end
            
        end


    xl=x{1,i};
    
       for l=1:Nk
        dist(l,k)=norm(x{1,i}(:,l)-x{1,i}(:,k));
        if (k==5) %disabled node
            DistanceDisabled(i,l)=norm(x{1,i}(1,l)-x{1,i}(1,k));
        end
        if (k==7)%general public
            DistanceK(i,l)=norm(x{1,i}(1,l)-x{1,i}(1,k));
        end
       end
       vind=find(dist(:,k)<R & dist(:,k)>0 );  
       %why W and xlk??
       W=[v{1,i}(1,k)/norm(v{1,i}(:,k)) -v{1,i}(2,k)/norm(v{1,i}(:,k)); v{1,i}(2,k)/norm(v{1,i}(:,k)) v{1,i}(1,k)/norm(v{1,i}(:,k)) ];
       p=1;    
       xlk=zeros(2,Nk);
       for l=vind'
            xlk(:,l)=W'*(xl(:,l)-x{1,i}(:,k)); 
            
            COFr{1,i}(:,k)=((COFvmax-v{1,i}(:,k))/(COFvmax-COFvmin))*(rmax-rmin)+rmin;  
            if (k==5) %double the distance for the disabled node
                COFr{1,i}(:,k)=COFr{1,i}(:,k)*2;
            end
            rnoise{1,i}(:,k)= COFr{1,i}(:,k); %add noise in the range of [-0.05,0.05]
            delta{1,i}(:,k)=delta{1,i}(:,k)+(xl(:,l)-x{1,i}(:,k))*(1-(rnoise{1,i}(1,k)/(norm(xl(:,l)-x{1,i}(:,k)))));
       end
       
       if xlk(1,:)==-abs(xlk(1,:))
           x{2,i}(1,k)=1;
       elseif xlk(2,:)==abs(xlk(2,:))
           x{2,i}(1,k)=2;
       elseif xlk(2,:)==-abs(xlk(2,:))
           x{2,i}(1,k)=3;
       else
           x{2,i}(1,k)=0;
       end
       
       if (va{1,i+1}(1,k)==0) && (va{1,i+1}(2,k)==0) 
           I=0;
       else
           I=1;
       end
       
     %%K1 faster than Kn
     maxDist1=1/(max(dist1(i,:)));
     minDist1=1/(min(dist1(i,:)));  
     Xdistance=1/dist1(i,k);
     COFdistance=((COFdistmax-COFdistmin)*(Xdistance-minDist1)/(maxDist1-minDist1))+COFdistmin;     
       
     delta{1,i}(:,k)=(1/(size(vind,1)-1))*delta{1,i}(:,k);
     vb{1,i+1}(:,k)=COFdistance*((1-lamda*I)*vg{1,i}(:,k))+gamma*delta{1,i}(:,k);
     
%      vb{1,i+1}(1,k)=COFdistance*((1-lamda*I)*vg{1,i}(1,k))+gamma*delta{1,i}(1,k);
%      vb{1,i+1}(2,k)=(1-lamda*I)*vg{1,i}(2,k)+gamma*delta{1,i}(2,k);

     v{1,i+1}(:,k)=lamda*I*(alpha*va{1,i+1}(:,k))+vb{1,i+1}(:,k);

     x{1,i+1}(:,k)=x{1,i}(:,k)+dt*v{1,i+1}(:,k);

     phi{1,i+1}(:,k)=((1-mu)*vg{1,i}(:,k)+mu*v{1,i+1}(:,k));
     
     
        if  (x{1,i+1}(1,k) >= wt2(1)-1)
            wt=wt3;

        
        elseif ((x{1,i+1}(1,k) >= wt1(1)-1) && (x{1,i+1}(1,k) < wt2(1)-1))
            wt=wt2;
        else
            wt=wt1;
        end
 

    end

    
    xl=x{1,i+1};
      for l=1:Nk
        dist(l,k)=norm(x{1,i+1}(:,l)-x{1,i+1}(:,k));
      end
       vind=find(dist(:,k)<R & dist(:,k)>0 );
 for k=1:Nk
  for l=vind'    
    vg{1,i+1}(:,k)=vg{1,i+1}(:,k)+phi{1,i+1}(:,l);
  end
  vg{1,i+1}(:,k)=(1/(size(vind,1)-1))*vg{1,i+1}(:,k);
 end

 
end

%% Graphical representation
%speed as scalar
for i=1:time
    for k=1:Nk
    ScalarSpeed(i,k)=sqrt(((v{1,i}(1,k))^2)+((v{1,i}(2,k))^2));
    end
end


for i=1:time
    %Display
    hold;
    Z=linspace(Xmin,Xmax, 100);
    scatter(x{1,i}(1,:),x{1,i}(2,:),'.','r'); %display of group of people in red
    if (k==5)
        scatter(x{1,i}(1,:),x{1,i}(2,:),'.','b');
    end

    scatter(wt(1,1),wt(2,1),'.','g'); %display of the target 
    scatter(Z,f1(Z),'.','b'); %display of circle 1
    scatter(Z,f2(Z),'.','b'); %display of circle 2

    axis([Xmin Xmax Ymin Ymax]);
    

    M(i) = getframe;
end

%%     %Display
    subplot(2,2,1);
   title("(a) i=0,Vaverage="+mean(ScalarSpeed(1,:))+"m/s")
    hold;
    Z=linspace(Xmin,Xmax, 200);
    scatter(x{1,1}(1,:),x{1,1}(2,:),'.','r'); %display of group of people in blue
    scatter(x{1,1}(1,5),x{1,1}(2,5),'.','k');
    scatter(wt(1,1),wt(2,1),'.','g'); %display of the target 
    scatter(Z,f1(Z),'.','b'); %display of circle 1
    scatter(Z,f2(Z),'.','b'); %display of circle 2
%      plot(MidPoints(1,:),MidPoints(2,:),'r--')
    axis([Xmin Xmax Ymin Ymax]);
%     hold off;
 
    subplot(2,2,2);
    title("(b) i=30,Vaverage="+mean(ScalarSpeed(30,:))+"m/s")
    hold;
    Z=linspace(Xmin,Xmax, 200);
    scatter(x{1,30}(1,:),x{1,30}(2,:),'.','r'); %display of group of people in blue
    scatter(x{1,30}(1,5),x{1,30}(2,5),'.','k');
    scatter(wt(1,1),wt(2,1),'.','g'); %display of the target 
    scatter(Z,f1(Z),'.','b'); %display of circle 1
    scatter(Z,f2(Z),'.','b'); %display of circle 2
%     plot(MidPoints(1,:),MidPoints(2,:),'r--')
    axis([Xmin Xmax Ymin Ymax]);

    subplot(2,2,3);
   title("(c) i=100,Vaverage="+mean(ScalarSpeed(100,:))+"m/s")
%     title ("(c) t=20")
    hold;
    Z=linspace(Xmin,Xmax, 200);
    scatter(x{1,100}(1,:),x{1,100}(2,:),'.','r'); %display of group of people in blue
    scatter(x{1,100}(1,5),x{1,100}(2,5),'.','k');
    scatter(wt(1,1),wt(2,1),'.','g'); %display of the target 
    scatter(Z,f1(Z),'.','b'); %display of circle 1
    scatter(Z,f2(Z),'.','b'); %display of circle 2
%     plot(MidPoints(1,:),MidPoints(2,:),'r--')
    axis([Xmin Xmax Ymin Ymax]);

    subplot(2,2,4);
   title("(d) i=200,Vaverage="+mean(ScalarSpeed(200,:))+"m/s")
%     title ("(d) t=60")
    hold;
    Z=linspace(Xmin,Xmax, 200);
    scatter(x{1,200}(1,:),x{1,200}(2,:),'.','r'); %display of group of people in blue
    scatter(x{1,200}(1,5),x{1,200}(2,5),'.','k');
    scatter(wt(1,1),wt(2,1),'.','g'); %display of the target 
    scatter(Z,f1(Z),'.','b'); %display of circle 1
    scatter(Z,f2(Z),'.','b'); %display of circle 2
%     plot(MidPoints(1,:),MidPoints(2,:),'r--')
    axis([Xmin Xmax Ymin Ymax]);
    
%% Gaphical representation of Euclidean Distance between each k node and target
for i=1:time
    for k=1:Nk    
        EuDist(i,k)=norm(wt-x{1,i}(:,k));
    end
end

for k=1:Nk
    hold on
    t=1:1:time;
%     plot(t,EuDist(:,k))
    if (k==5)
        plot(t,EuDist(:,k), 'r')
    else
        plot(t,EuDist(:,k), 'b')
    end
end
