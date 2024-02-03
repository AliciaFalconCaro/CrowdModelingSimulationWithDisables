%% Crowd Modeling Simulation 

% A smaller dt increases the accuracy of the system, but also increases the simulation time.

clear all;close all;clc;

% Window frame of the simulation
Xmin=-2;        
Ymin=-4;
Xmax=20;
Ymax=3;

Nk=40;       % no. of nodes
dt=0.05;     % Time increment / time step
time=250;    % Number of iterations      

% Speed Coefficients and adjustable parameters
lamda=0.5;  
beta=1;  
gamma=1;  
mu=0.5;  
alpha=0.5;

R=15;       % Neighborhood radius
rmax=2;     % Neighbors' distance from each other
rmin=1.5;
Lcmin=1;    % Max and Min chord
Lcmax=10;
Vmin=0.6;   % Max and Min speed
Vmax=3.6; 
%Vo=2.1;     % Average speed
Cki=0;      % C_{k,i}=v_{k,i}^{c}

% Distance between k and target, used to make some nodes go faster than others (tail effect)
TailDistanceMin=0.1; 
TailDistanceMax=4; 

va=cell(1,time+1);
vb=cell(1,time+1);
delta=cell(1,time);
v=cell(1,time+1);
v{1,1}=zeros(2,Nk);
phi=cell(1,time+1);
vg=cell(1,time+1);
vg{1,1}=zeros(2,Nk);
Lck=zeros(time, Nk);
rnoise=cell(1,time);
rk=cell(1,time);
w{1,1}=randn(2,Nk);     % Random initial location of nodes

% Target(s):
t1=[6;1];
t2=[13;-1];
t3=[18; -3];
t=t1;

COFv=0; 

% Borders 
syms Z1 Z2
f11= 2*cos(((2*pi)/30)*(Z1-1.5))+1;
f22= 4*sin(((2*pi)/40)*(Z2+4))-3.5;

f1= @(Z)2*cos(((2*pi)/30)*(Z-1.5))+1;
f2= @(Z)4*sin(((2*pi)/40)*(Z+4))-3.5;


% To calculate the chords, we obtain all the possible centres of the circles
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

% To ensure all the nodes are randomly located within the walls at the start of the simulation
for i=1:Nk
    xCoord=w{1,1}(1,i);
    yCoord=w{1,1}(2,i);
    if (yCoord >= f1(xCoord))
        w{1,1}(2,i)=f1(xCoord)-1;
    end
    if (yCoord <= f2(xCoord))
        w{1,1}(2,i)=f2(xCoord)+1;
    end
end

%% Main program
for i=1:time
    vg{1,i+1}=zeros(2,Nk);
    delta{1,i}=zeros(2,Nk);
    rnoise{1,i}=zeros(2,Nk);
    rk{1,i}=zeros(2,Nk);

    for k=1:Nk    
     DistanceToTarget(i,k)=norm(t-w{1,i}(:,k));
    end
           
    for k=1:Nk
        % Calculate chord through node k at time i
        for j=1:length(MidPoints)
            MidPoint_j=MidPoints(:,j);
            if (MidPoint_j(1)==round(w{1,i}(1,k),1))
                PointC=MidPoints(:,j);
                PointD=MidPoints(:,j-1);
                coefficients = polyfit(PointC, PointD, 1); %eq between Point C and D
                fPerpendicularMidPointAB1=(-1/coefficients(1))*(Z1-w{1,i}(1,k))+w{1,i}(2,k);
                fPerpendicularMidPointAB2=(-1/coefficients(1))*(Z2-w{1,i}(1,k))+w{1,i}(2,k);
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
        
        Cki(i,k)=((Lcmax - Lck(i,k))/(Lcmax-Lcmin))*(Vmax-Vmin)+Vmin;

        if (k==5) % Half the speed for the disabled node
            Cki(i,k)=Cki(i,k)/2;
        end

         XCoordinate_w=w{1,i}(1,k);
         pki_UpperWall=[XCoordinate_w;(f1(XCoordinate_w)-alpha)];
         pki_LoweWall=[XCoordinate_w;(f2(XCoordinate_w)+alpha)];
 
        if (w{1,i}(2,k)<(f1(XCoordinate_w)) && w{1,i}(2,k)>(f2(XCoordinate_w))) % node k within walls
            va{1,i+1}(:,k)=Cki(i,k)*((t-w{1,i}(:,k))./norm(t-w{1,i}(:,k)));   
        else 
            DistanceToUpperWall = norm(pki_UpperWall-w{1,i}(:,k));
            if DistanceToUpperWall>3
               DistanceToUpperWall=3; 
            end
            if DistanceToUpperWall<1
               DistanceToUpperWall=1; 
            end
            DistanceToLowerWall = norm(pki_LoweWall-w{1,i}(:,k));
            if DistanceToLowerWall>3
               DistanceToLowerWall=3; 
            end
            if DistanceToLowerWall<1
               DistanceToLowerWall=1; 
            end

            if (w{1,i}(2,k)>=(f1(XCoordinate_w))) 
               va{1,i+1}(:,k)=-Cki(i,k)*((R/DistanceToUpperWall)-1)*(w{1,i}(:,k)-pki_UpperWall);
            elseif (w{1,i}(2,k)<=(f2(XCoordinate_w)))
               va{1,i+1}(:,k)=-Cki(i,k)*((R/DistanceToLowerWall)-1)*(w{1,i}(:,k)-pki_LoweWall);
            end 
        end

        wl=w{1,i};    
        for l=1:Nk
            distance_kl(l,k)=norm(w{1,i}(:,l)-w{1,i}(:,k));
            if (k==5) % Disabled node
                DistanceDisabled(i,l)=norm(w{1,i}(1,l)-w{1,i}(1,k));
            end
            if (k==7)% General public
                DistanceK(i,l)=norm(w{1,i}(1,l)-w{1,i}(1,k));
            end
        end
        vind=find(distance_kl(:,k)<R & distance_kl(:,k)>0 );  
        W=[v{1,i}(1,k)/norm(v{1,i}(:,k)) -v{1,i}(2,k)/norm(v{1,i}(:,k)); v{1,i}(2,k)/norm(v{1,i}(:,k)) v{1,i}(1,k)/norm(v{1,i}(:,k)) ];
        %p=1;    
        wlk=zeros(2,Nk);
    
        for l=vind'
            wlk(:,l)=W'*(wl(:,l)-w{1,i}(:,k));  
            rk{1,i}(:,k)=((Vmax-v{1,i}(:,k))/(Vmax-Vmin))*(rmax-rmin)+rmin;  
        
            if (k==5) % Double the distance for the disabled node
                rk{1,i}(:,k)=rk{1,i}(:,k)*2;
            end
        
            rnoise{1,i}(:,k)= rk{1,i}(:,k)+(rand()-0.05); % Add noise here if necessary
            delta{1,i}(:,k)=delta{1,i}(:,k)+(wl(:,l)-w{1,i}(:,k))*(1-(rnoise{1,i}(1,k)/(norm(wl(:,l)-w{1,i}(:,k)))));
        end
       
        if wlk(1,:)==-abs(wlk(1,:))
            w{2,i}(1,k)=1;
        elseif wlk(2,:)==abs(wlk(2,:))
            w{2,i}(1,k)=2;
        elseif wlk(2,:)==-abs(wlk(2,:))
            w{2,i}(1,k)=3;
        else
            w{2,i}(1,k)=0;
        end
       
        if (va{1,i+1}(1,k)==0) && (va{1,i+1}(2,k)==0) 
            I=0;
        else
            I=1;
        end
       
        % Tail effect
        MaxNodesDist=1/(max(DistanceToTarget(i,:)));
        MinNodesDist=1/(min(DistanceToTarget(i,:)));  
        TailDistance=((TailDistanceMin-TailDistanceMax)*((1/DistanceToTarget(i,k))-MinNodesDist)/(MaxNodesDist-MinNodesDist))+TailDistanceMax;      
     
        delta{1,i}(:,k)=(1/(size(vind,1)-1))*delta{1,i}(:,k);
        vb{1,i+1}(:,k)=TailDistance*((1-lamda*I)*vg{1,i}(:,k))+gamma*delta{1,i}(:,k);
        v{1,i+1}(:,k)=lamda*I*(beta*va{1,i+1}(:,k))+vb{1,i+1}(:,k);
        w{1,i+1}(:,k)=w{1,i}(:,k)+dt*v{1,i+1}(:,k);
        phi{1,i+1}(:,k)=((1-mu)*vg{1,i}(:,k)+mu*v{1,i+1}(:,k));
     
        if  (w{1,i+1}(1,k) >= t2(1)-1)
            t=t3;
        elseif ((w{1,i+1}(1,k) >= t1(1)-1) && (w{1,i+1}(1,k) < t2(1)-1))
            t=t2;
        else
            t=t1;
        end
    end

    % Combine step
    wl=w{1,i+1};
    for l=1:Nk
        distance_kl(l,k)=norm(w{1,i+1}(:,l)-w{1,i+1}(:,k));
    end
    vind=find(distance_kl(:,k)<R & distance_kl(:,k)>0 );
    for k=1:Nk
        for l=vind'    
            vg{1,i+1}(:,k)=vg{1,i+1}(:,k)+phi{1,i+1}(:,l);
        end
        vg{1,i+1}(:,k)=(1/(size(vind,1)-1))*vg{1,i+1}(:,k);
    end
end
