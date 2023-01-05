clc
clear
tic

% matlabpool('open',2);

% load PORE_NETWORK_DATA#H2cm;
load Berea_S4;
[N_PORE_DATA,N_PORE_DATA2]=size(PORE);
[N_THROAT_DATA,N_THROAT_DATA2]=size(THROAT);

%% Defining variables
PORE_STATUS=zeros(N_PORE_DATA,1);
THROAT_STATUS=zeros(N_THROAT_DATA,1);
PORE_CONDUCTANCE=zeros(N_PORE_DATA,1);
THROAT_CONDUCTANCE=zeros(N_THROAT_DATA,1);
THROAT_CONDUCTANCE2=zeros(N_THROAT_DATA,1);
AVERAGE_THROAT_CONDUCTANCE=zeros(N_THROAT_DATA,1);

A=zeros(N_PORE_DATA,N_PORE_DATA);
C=zeros(N_PORE_DATA,1);

PORE_INLET_STATUS=zeros(N_PORE_DATA,1);
PORE_OUTLET_STATUS=zeros(N_PORE_DATA,1);
THROAT_INLET_STATUS=zeros(N_THROAT_DATA,1);
THROAT_OUTLET_STATUS=zeros(N_THROAT_DATA,1);

THROAT_Q=zeros(N_THROAT_DATA,1);
QIJ=zeros(N_PORE_DATA,100);

PORE_AREA=zeros(N_PORE_DATA,1);
THROAT_AREA=zeros(N_THROAT_DATA,1);
THROAT_AVERAGE_AREA=zeros(N_THROAT_DATA,1);

MI=zeros(N_PORE_DATA,1);
NIJ=zeros(N_THROAT_DATA,1);
BIJ=zeros(N_THROAT_DATA,1);




%% Boundary condition
Source_Pressure=30; %pa
Outlet_Pressure=1; %pa

WaterViscosity=8.90e-4; %pa.s (SI)

%% Initialization
for i=1:N_PORE_DATA                    % This for loop determines what is the shape of the pore (circle, square, or triangle)
    if abs(PORE(i,7)-0.079577)<0.001
        PORE_STATUS(i)=1;               %1 for circle
    elseif abs(PORE(i,7)-0.0625)<0.001
        PORE_STATUS(i)=0;               %0 for square
    else
        PORE_STATUS(i)=-1;              %-1 for triangle
    end
end





for i=1:N_PORE_DATA
    if PORE_STATUS(i)==1
	    PORE_CONDUCTANCE(i)=0.5*(PORE(i,6)^4)/16/WaterViscosity/PORE(i,7); 
	
    elseif PORE_STATUS(i)==0
		PORE_CONDUCTANCE(i)=0.5623*(PORE(i,6)^4)/16/WaterViscosity/PORE(i,7);  
	
    else
		PORE_CONDUCTANCE(i)=0.6*(PORE(i,6)^4)/16/WaterViscosity/PORE(i,7);  
	
    end
end


% WRONG
% for i=1:N_PORE_DATA
%     if PORE_STATUS(i)==1
% 	    PORE_CONDUCTANCE(i)=0.5*(PORE(i,6)^4)*pi/4/WaterViscosity;  %Pore_radius^4 *pi/8
% 	
%     elseif PORE_STATUS(i)==0
% 		PORE_CONDUCTANCE(i)=0.5623*(PORE(i,6)^4)*pi/4/WaterViscosity;  %0.5623*Pore_radius^4 *pi/4
% 	
%     else
% 		PORE_CONDUCTANCE(i)=0.6*(PORE(i,6)^4)*pi/4/WaterViscosity;  %0.6*Pore_radius^4 *pi/4
% 	
%     end
% end



for i=1:N_THROAT_DATA                    % This for loop determines what is the shape of the throat (circle, square, or triangle)
    if abs(THROAT(i,4)-0.079577)<0.001
        THROAT_STATUS(i)=1;               %1 for circle
    elseif abs(THROAT(i,4)-0.0625)<0.001
        THROAT_STATUS(i)=0;               %0 for square
    else
        THROAT_STATUS(i)=-1;              %-1 for triangle
    end
end

%for i=1:N_THROAT_DATA
%    THROAT_CONDUCTANCE(i)=(THROAT(i,3)^4)*pi/8; %Throat_radius^4 *pi/8
%end

for i=1:N_THROAT_DATA
    if THROAT_STATUS(i)==1
	    THROAT_CONDUCTANCE(i)=0.5*(THROAT(i,3)^4)/16/WaterViscosity/THROAT(i,4);  %Conductance of a circle Throat_radius^4 *pi/8
	elseif THROAT_STATUS(i)==0
		THROAT_CONDUCTANCE(i)=0.5623*(THROAT(i,3)^4)/16/WaterViscosity/THROAT(i,4);  %conductance of square: 0.5623*Throat_radius^4 *pi/4
    else
		THROAT_CONDUCTANCE(i)=0.6*(THROAT(i,3)^4)/16/WaterViscosity/THROAT(i,4);  %Conductance of triangle 0.6*Throat_radius^4 *pi/4
	
    end
end




% WRONG
% for i=1:N_THROAT_DATA
%     if THROAT_STATUS(i)==1
% 	    THROAT_CONDUCTANCE(i)=0.5*(THROAT(i,3)^4)*pi/4/WaterViscosity;  %Conductance of a circle Throat_radius^4 *pi/8
% 	elseif THROAT_STATUS(i)==0
% 		THROAT_CONDUCTANCE(i)=0.5623*(THROAT(i,3)^4)*pi/4/WaterViscosity;  %conductance of square: 0.5623*Throat_radius^4 *pi/4
% 	else
% 		THROAT_CONDUCTANCE(i)=0.6*(THROAT(i,3)^4)*pi/4/WaterViscosity;  %Conductance of triangle 0.6*Throat_radius^4 *pi/4
% 	
%     end
% end


for i=1:N_PORE_DATA                       % This for loop determines if the pore is connected to the inlet or outlet or none.
    if CONNECTIONDATA(i,PORE(i,4)+1)==1
        PORE_INLET_STATUS(i)=1;
    else
        PORE_INLET_STATUS(i)=0;
    end
    
    if CONNECTIONDATA(i,PORE(i,4)+2)==1
        PORE_OUTLET_STATUS(i)=1;
    else
        PORE_OUTLET_STATUS(i)=0;
    end
end

for i=1:N_THROAT_DATA
    if (THROAT(i,1)==-1) || (THROAT(i,2)==-1)
        THROAT_INLET_STATUS(i)=1;
    else
        THROAT_INLET_STATUS(i)=0;
    end
    
    if (THROAT(i,1)==0) || (THROAT(i,2)==0)
        THROAT_OUTLET_STATUS(i)=1;
    else
        THROAT_OUTLET_STATUS(i)=0;
    end    
    
end




for i=1:N_THROAT_DATA               % This for loop calculates the average conductance of a throat using conductance of connected pores
    if (THROAT(i,1)==-1) || THROAT(i,1)==0
        AVERAGE_THROAT_CONDUCTANCE(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
            ( (THROAT(i,7)/(PORE_CONDUCTANCE(THROAT(i,2)))) + ...
            (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) );  %Here it assumed that inlet or outlet imaginary nodes have infinite conductance. So the ratio of pore length devided by pore conductance will be zero
        
        THROATEFFECTIVELENGTH(i,1)=THROAT(i,6)+THROAT(i,7)+THROAT(i,8);
        GavPerLtotal(i)=1/( (THROAT(i,7)/(PORE_CONDUCTANCE(THROAT(i,2)))) + (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) ); 
        
        
    elseif (THROAT(i,2)==-1) || THROAT(i,2)==0
        AVERAGE_THROAT_CONDUCTANCE(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
            ( (THROAT(i,6)/(PORE_CONDUCTANCE(THROAT(i,1)))) + ...
            (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) );
        
        THROATEFFECTIVELENGTH(i,1)=THROAT(i,6)+THROAT(i,7)+THROAT(i,8);
        GavPerLtotal(i)=1/( (THROAT(i,6)/(PORE_CONDUCTANCE(THROAT(i,1))))  + (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) ); 
        
        
    else
        
        AVERAGE_THROAT_CONDUCTANCE(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
            ( (THROAT(i,6)/(PORE_CONDUCTANCE(THROAT(i,1)))) + ...
             (THROAT(i,7)/(PORE_CONDUCTANCE(THROAT(i,2)))) + ...
            (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) );
        THROATEFFECTIVELENGTH(i,1)=THROAT(i,6)+THROAT(i,7)+THROAT(i,8);
        GavPerLtotal(i)=1/( (THROAT(i,6)/(PORE_CONDUCTANCE(THROAT(i,1))))  + (THROAT(i,7)/(PORE_CONDUCTANCE(THROAT(i,2)))) + (THROAT(i,8)/(THROAT_CONDUCTANCE(i))) ); 
        
        
    end
end



%% determining the matrix A and C in A*P=C

for i=1:N_PORE_DATA
    
    if PORE_INLET_STATUS(i)==1
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)~=-1
                A(i,CONNECTIONDATA(i,a))=-AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
            end
        end
        for a=1:PORE(i,4)
            A(i,i)=A(i,i)+ AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
        end
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)==-1
                C(i)=Source_Pressure * AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
            end
        end        
        
        
    elseif PORE_OUTLET_STATUS(i)==1
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)~=0
                A(i,CONNECTIONDATA(i,a))=-AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
            end
        end
        for a=1:PORE(i,4)
            A(i,i)=A(i,i)+ AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
        end
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)==0
                C(i)=Outlet_Pressure * AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
            end
        end 
        
    else
        for a=1:PORE(i,4)
            A(i,CONNECTIONDATA(i,a))=-AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
        end
        for a=1:PORE(i,4)
            A(i,i)=A(i,i)+ AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5);
        end
        C(i)=0;

    end
end




%% Solving the system of equations to obtain pressure of each node

% P=inv(A)*C;

[ACorrected1,CCorrected1,P]=PressureSolve(A,C);

%% Determining flow through each throat

for i=1:N_THROAT_DATA
    if THROAT(i,1) ==-1
        P_Upstream=Source_Pressure;
        P_Downstream=P(THROAT(i,2));
    elseif THROAT(i,1)==0
        P_Upstream=Outlet_Pressure;
        P_Downstream=P(THROAT(i,2));
    elseif THROAT(i,2)==-1
        P_Upstream=P(THROAT(i,1));
        P_Downstream=Source_Pressure;  
    elseif THROAT(i,2)==0
        P_Upstream=P(THROAT(i,1));
        P_Downstream=Outlet_Pressure;
    else
        P_Upstream=P(THROAT(i,1));
        P_Downstream=P(THROAT(i,2));
    end    
    
    THROAT_Q(i)=(AVERAGE_THROAT_CONDUCTANCE(i)/THROAT(i,5))*(P_Upstream-P_Downstream);
    
    
    Qin(i)=abs(THROAT_INLET_STATUS(i)*THROAT_Q(i));
    Qout(i)=abs(THROAT_OUTLET_STATUS(i)*THROAT_Q(i));
    
end

Q_INPUT=sum(Qin);
Q_OUTPUT=sum(Qout);

Material_Balance=Q_OUTPUT-Q_INPUT

% scatter3(PORE(:,1),PORE(:,2),PORE(:,3), 300000*PORE(:,6),P,'filled')

GuessTime=3*(sum(THROAT(:,9))+sum(PORE(:,5)))/Q_OUTPUT
    
%% To calculate Qij in each node

% Cross-check after calculation of QIJ is to sum each row. The elements of
% each row are inflowing and outflowing streams of a node and their sum is
% zero.

for i=1:N_PORE_DATA
    
    if PORE_INLET_STATUS(i)==1
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)==-1
                QIJ(i,a)=(AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5))*(P(i)-Source_Pressure);
            else
                QIJ(i,a)=(AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5))*(P(i)-P(CONNECTIONDATA(i,a)));
            end
        end
 
    elseif PORE_OUTLET_STATUS(i)==1
        for a=1:PORE(i,4)
            if CONNECTIONDATA(i,a)==0
                QIJ(i,a)=(AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5))*(P(i)-Outlet_Pressure);
            else
                QIJ(i,a)=(AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5))*(P(i)-P(CONNECTIONDATA(i,a)));
            end
        end
        
        
    else
        for a=1:PORE(i,4)
            QIJ(i,a)=(AVERAGE_THROAT_CONDUCTANCE(CONNECTIONDATA(i,PORE(i,4)+2+a))/THROAT(CONNECTIONDATA(i,PORE(i,4)+2+a),5))*(P(i)-P(CONNECTIONDATA(i,a)));
        end
    end
    
end

% for i=1:N_PORE_DATA
%     
%     if PORE_INLET_STATUS(i)==1
%         for a=1:PORE(i,4)
%             if CONNECTIONDATA(i,a)==-1
%                 QIJ(i,a)=(GavPerLtotal(CONNECTIONDATA(i,PORE(i,4)+2+a)))*(P(i)-Source_Pressure);
%             else
%                 QIJ(i,a)=(GavPerLtotal(CONNECTIONDATA(i,PORE(i,4)+2+a)))*(P(i)-P(CONNECTIONDATA(i,a)));
%             end
%         end
%  
%     elseif PORE_OUTLET_STATUS(i)==1
%         for a=1:PORE(i,4)
%             if CONNECTIONDATA(i,a)==0
%                 QIJ(i,a)=(GavPerLtotal(CONNECTIONDATA(i,PORE(i,4)+2+a)))*(P(i)-Outlet_Pressure);
%             else
%                 QIJ(i,a)=(GavPerLtotal(CONNECTIONDATA(i,PORE(i,4)+2+a)))*(P(i)-P(CONNECTIONDATA(i,a)));
%             end
%         end
%         
%         
%     else
%         for a=1:PORE(i,4)
%             QIJ(i,a)=(GavPerLtotal(CONNECTIONDATA(i,PORE(i,4)+2+a)))*(P(i)-P(CONNECTIONDATA(i,a)));
%         end
%     end
%     
% end

%% These section is to calculate cross section area of pores and throats

% 
for i=1:N_PORE_DATA                    
    PORE_AREA(i)=(PORE(i,6))^2/(4*PORE(i,7)); % r^2/(4G)
end

for i=1:N_THROAT_DATA                    
    THROAT_AREA(i)=(THROAT(i,3))^2/(4*THROAT(i,4)); % r^2/(4G)
end

for i=1:N_THROAT_DATA
    if THROAT(i,1) == -1
        THROAT_AVERAGE_AREA(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
                           (  (THROAT(i,6)/1.5e-5 ) ...
                             +(THROAT(i,7)/PORE_AREA(THROAT(i,2)) ) ...
                             +(THROAT(i,8)/THROAT_AREA(i)) );
    elseif THROAT(i,2) == 0
        THROAT_AVERAGE_AREA(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
                           (  (THROAT(i,6)/PORE_AREA(THROAT(i,1)) ) ...
                             +(THROAT(i,7)/1.5e-5 ) ...
                             +(THROAT(i,8)/THROAT_AREA(i)) );
    else
        THROAT_AVERAGE_AREA(i)=(THROAT(i,6)+THROAT(i,7)+THROAT(i,8))/...
                           (  (THROAT(i,6)/PORE_AREA(THROAT(i,1)) ) ...
                             +(THROAT(i,7)/PORE_AREA(THROAT(i,2)) ) ...
                             +(THROAT(i,8)/THROAT_AREA(i)) );
    end
end



%% Parameter definition for Solute transport

deltaT=0.001;                    %second
Time=GuessTime;                        %second
NOTS=ceil(Time/deltaT);


% D=1e-9; %Diffusion coefficient; %m2/sec
% D=0;
for i=1:N_THROAT_DATA
    if  THROAT_INLET_STATUS(i)==1
        D(i)=0;
    elseif THROAT_OUTLET_STATUS(i)==1 
        D(i)=0;
    else
        D(i)=1e-9;
    end    

              
  
end



A_Solute=zeros(N_PORE_DATA,N_PORE_DATA);
Z_Solute=zeros(N_PORE_DATA,1);


CONC_P=zeros(N_PORE_DATA,1);
CONC_T=zeros(N_THROAT_DATA,1);



%% Initial and boundary conditions for concentration

CONC_P(:,:)=0;
CONC_T(:,:)=0;

Source_Concentration=1;
OUTLET_Concentration=0;

%% Initialization


for i=1:N_PORE_DATA                    
    MI(i)=deltaT/(PORE(i,5)); % dt/Volume
end

for i=1:N_THROAT_DATA
    if THROAT(i,9)==0 
        NIJ(i)=deltaT/(1e-12); % dt/Volume
    else
        NIJ(i)=deltaT/(THROAT(i,9)); % dt/Volume
    end
end


for i=1:N_THROAT_DATA
    BIJ(i)=2*D(i)*THROAT_AVERAGE_AREA(i)/(THROAT(i,8)); % 2*D*A/L
end



%% SIMULATION
for timestep=1:NOTS
A_Solute(:,:)=0;
Z_Solute(:,:)=0; %Resetting the matrices of the equation system
    
    for i=1:N_PORE_DATA
        A_Solute(i,i)=1;
        if PORE_INLET_STATUS(i)==1
            for a=1:PORE(i,4)
                if CONNECTIONDATA(i,a)~=-1
                    A_Solute(i,CONNECTIONDATA(i,a))= -MI(i)*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0))^2 ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
                end
            end
            for a=1:PORE(i,4)
                A_Solute(i,i)=A_Solute(i,i)+ MI(i)*max(QIJ(i,a),0) + MI(i)*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    - MI(i)* NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*(abs(QIJ(i,a))+BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)))...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
                
                Z_Solute(i)= Z_Solute(i)+ ...
                    CONC_T(CONNECTIONDATA(i,PORE(i,4)+2+a)) * ...
                    MI(i) * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0)) ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
            end
            Z_Solute(i)=Z_Solute(i)+CONC_P(i) ...
                + Source_Concentration * MI(i)*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0))^2  ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
            
            
        elseif PORE_OUTLET_STATUS(i)==1
            for a=1:PORE(i,4)
                if CONNECTIONDATA(i,a)~=0
                    A_Solute(i,CONNECTIONDATA(i,a))= -MI(i)*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0))^2 ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
                end
            end
            for a=1:PORE(i,4)
                A_Solute(i,i)=A_Solute(i,i)+ MI(i)*max(QIJ(i,a),0) + MI(i)*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    - MI(i)* NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*(abs(QIJ(i,a))+BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)))...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
                
                Z_Solute(i)= Z_Solute(i)+ ...
                    CONC_T(CONNECTIONDATA(i,PORE(i,4)+2+a)) * ...
                    MI(i) * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0)) ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
            end
            Z_Solute(i)=Z_Solute(i) +  CONC_P(i) ...
                + OUTLET_Concentration * MI(i)*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                   * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0))^2  ...
                   /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
              
        else

            for a=1:PORE(i,4)
                A_Solute(i,CONNECTIONDATA(i,a))= -MI(i)*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0))^2 ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));

                A_Solute(i,i)=A_Solute(i,i)+ MI(i)*max(QIJ(i,a),0) + MI(i)*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)) ...
                    - MI(i)* NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*(abs(QIJ(i,a))+BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)))...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
                
                Z_Solute(i)= Z_Solute(i)+ ...
                    CONC_T(CONNECTIONDATA(i,PORE(i,4)+2+a)) * ...
                    MI(i) * (BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))-min(QIJ(i,a),0)) ...
                    /(1 + NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*abs(QIJ(i,a)) + 2*NIJ(CONNECTIONDATA(i,PORE(i,4)+2+a))*BIJ(CONNECTIONDATA(i,PORE(i,4)+2+a)));
            end
            Z_Solute(i)=Z_Solute(i) +  CONC_P(i);
                   
        end
        
    end
    
    % Solving by inversing matric A_Solute and multiplying it to Z_Solute.

    CONC_P=A_Solute\Z_Solute;
    
    
%     for i=1:N_PORE_DATA
%         if CONC_P(i)>1
%             CONC_P(i)=1;
%         end
%     end

    
%     CONC_PORE(:,timestep)=CONC_P;
    
    % Solving the explicit equation for throats
    
    for i=1:N_THROAT_DATA
        
        if THROAT(i,1) ==-1
            CONC_T(i)= ( ...
            Source_Concentration *  ( NIJ(i) *max(THROAT_Q(i),0)  + NIJ (i)*BIJ(i))...
            -CONC_P(THROAT(i,2)) *  ( NIJ(i) *min(THROAT_Q(i),0)  - NIJ (i)*BIJ(i))...
            + CONC_T(i))...
            /(1 + NIJ(i)*abs(THROAT_Q(i)) + 2*NIJ(i)*BIJ(i));
        
        elseif THROAT(i,1)==0
        
            CONC_T(i)= ( ...
            OUTLET_Concentration *  ( NIJ(i) *max(THROAT_Q(i),0)  + NIJ (i)*BIJ(i))...
            -CONC_P(THROAT(i,2)) *  ( NIJ(i) *min(THROAT_Q(i),0)  - NIJ (i)*BIJ(i))...
            + CONC_T(i))...
            /(1 + NIJ(i)*abs(THROAT_Q(i)) + 2*NIJ(i)*BIJ(i));
       
        elseif THROAT(i,2)==-1
            
            CONC_T(i)= ( ...
            CONC_P(THROAT(i,1)) *  ( NIJ(i) *max(THROAT_Q(i),0)  + NIJ (i)*BIJ(i))...
            - Source_Concentration *  ( NIJ(i) *min(THROAT_Q(i),0)  - NIJ (i)*BIJ(i))...
            + CONC_T(i))...
            /(1 + NIJ(i)*abs(THROAT_Q(i)) + 2*NIJ(i)*BIJ(i));
        
    
        elseif THROAT(i,2)==0
            
            CONC_T(i)= ( ...
             CONC_P(THROAT(i,1))*  ( NIJ(i) *max(THROAT_Q(i),0)  + NIJ (i)*BIJ(i))...
            -OUTLET_Concentration *  ( NIJ(i) *min(THROAT_Q(i),0)  - NIJ (i)*BIJ(i))...
            + CONC_T(i))...
            /(1 + NIJ(i)*abs(THROAT_Q(i)) + 2*NIJ(i)*BIJ(i));
        
        else
        CONC_T(i)= ( ...
            CONC_P(THROAT(i,1)) *  ( NIJ(i) *max(THROAT_Q(i),0)  + NIJ (i)*BIJ(i))...
            -CONC_P(THROAT(i,2)) *  ( NIJ(i) *min(THROAT_Q(i),0)  - NIJ (i)*BIJ(i))...
            + CONC_T(i))...
            /(1 + NIJ(i)*abs(THROAT_Q(i)) + 2*NIJ(i)*BIJ(i));
        end
    end
    
    
    for i=1:N_THROAT_DATA
        Solute_Mass_OUT(i)=abs(THROAT_OUTLET_STATUS(i)*THROAT_Q(i)*CONC_T(i));
    end
    
    OUTLET_Concentration=sum(Solute_Mass_OUT)/Q_OUTPUT;
%     CONC_THROAT(:,timestep)=CONC_T;
    OUTLET_CONC(timestep,1)=OUTLET_Concentration;
    TIME(timestep)=timestep*deltaT;
%     timestep

    if mod(timestep, 2000)==1
        h=figure('Visible','off');
        
        scatter3(PORE(:,1),PORE(:,2),PORE(:,3), 3000000*PORE(:,6),CONC_P,'filled')
        axis equal
%            set(gca,'Zlim',[-0.01 0.01]);
%           set(gca,'Ylim',[-0.01 0.01]);
        colorbar;
        
        NameOfOutputPNGFile=strcat('cocentration',num2str(timestep),'.png');
        NameOfsavefile=strcat('saving time=',num2str(timestep),'.mat');
        saveas(h,NameOfOutputPNGFile,'png');
        save(NameOfsavefile);
        
        
    end
    
    
end


CONC_P
CONC_T
% Q_INPUT
Q_OUTPUT


toc


n=1;
for i=1:1000000
    if mod(i,100)==1
        x(n,1)=TIME(i);
        y(n,1)=OUTLET_CONC(i);
        n=n+1;
    end
end




% % matlabpool('close');
% 
% % 
% for i=1:N_PORE_DATA
%     for j=1:100
%         QIJP(i,j)=max(QIJ(i,j),0);
%     end
%     QINP(i,1)=sum(QIJP(i,:));
%     STABLEPORE(i,1)=PORE(i,5)/(QINP(i)); % dt/Volume
% end
% 
% for i=1:N_THROAT_DATA
%     if THROAT(i,9)==0
%         THROAT_VOLUME(i)=1e-12;
%     else
%         THROAT_VOLUME(i)=THROAT(i,9);
%     end
%     STABLETHROAT(i,1)=THROAT_VOLUME(i)/(abs(THROAT_Q(i))); % dt/Volume
% end
% 
% 
% 
% % 
for i=1:N_PORE_DATA
    MASSPORE(i)=PORE(i,5)*CONC_P(i);
end

for i=1:N_THROAT_DATA
    if THROAT(i,9)==0
        THROAT_VOLUME(i)=1e-11;
    else
        THROAT_VOLUME(i)=THROAT(i,9);
    end
    MASSTHROAT(i)=THROAT_VOLUME(i)*CONC_T(i);
end

SOLUTEMASS=(sum(MASSTHROAT(:))+sum(MASSPORE(:)))


for i=1:N_PORE_DATA
    sumQ(i,1)=sum(QIJ(i,:));
end




    