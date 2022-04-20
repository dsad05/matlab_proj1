function cz3_integr
    clear;
    close all;
    clc;
    syms t;

    T=3;
    F1=5;
    F2=50;

    syms t;
    E=@(t) 1;

    %{ 
    E=@(t) 1;
    E=@(t) period(t,T);
    E=@(t) 240*sin(t);
    E=@(t) 210*sin(2*pi*F1*t);
    E=@(t) 120*sin(2*pi*F2*t);
    %}

    %%Select step
    %h=0.0005;
    h=0.005;

    function e = period(t,T)
        p=rem(t,T);
        if(p>= T/2)
            e=0;
        else
            e=120;
        end
    end

    a = 0;
    b = 30;
    t = a:h:b;      %[s]
    y=[0;           %[A]
       0;           %[A]   
       0];          %[V]
    R1 = 0.1;       %[Ohm]
    R2 = 10;        %[Ohm]
    C = 0.5;        %[F]
    L1 = 3;         %[H]
    L2 = 5;         %[H]
    
    i1_all=zeros(1,length(t));
    i2_all=zeros(1,length(t));
    uC_all=zeros(1,length(t));
    P_inst=zeros(1,length(t));

        y=modified_euler(t, y, h);
        for k=1:length(t)
            i1_all(1,:)=y(1,:);
            i2_all(1,:)=y(2,:);
            uC_all(1,:)=y(3,:);
            P_inst(1,:)=(y(1,:).^2).*R1+(y(2,:).^2).*R2; %instantaneous power
        end

    I_r = int_r(P_inst, h);
    I_r_sum=sum(I_r)
    I_p = int_p(P_inst, h);
    I_p_sum=sum(I_p)
    P_vecsum=zeros(1,length(t));
    P_vecsum(1,1)=I_p(1);
    for ii=2:length(t)
        P_vecsum(ii)=I_p(ii)+P_vecsum(ii-1);
    end


   fig=figure;  
   fig.WindowState = 'maximized';
    subplot(3,2,2);                    
        plot(t, P_inst(1,:))
        grid on;
        legend('Instantaneous power');
        ylabel('P[W]');
        xlabel('t[s]');
        title('P(t)');
% hold on;

%         subplot(3,2,4);                    
%         plot(t,P_vecsum)
%         grid on;
%         legend('sum of instantaneous power');
%         ylabel('P[W]');
%         xlabel('t[s]');
%         title('P(t)');

    subplot(3,2,3.5);              
        plot(t, i1_all(1,:), ...
            t, i2_all(1,:),'--')
        grid on;
        legend('I1[A]', ...
            'I2[A]');
        ylabel('I[A]');
        xlabel('t[s]');
        title('I(t)');
        
    subplot(3,2,5);                   
        plot(t, uC_all(1,:))
        grid on;
        legend('UC[V]');
        ylabel('UC[V]');
        xlabel('t[s]');
        title('UC(t)');

    function dy = f(t,y)
        e=E(t);
        i1=y(1);
        i2=y(2);
        uC=y(3);
        uL1 =   [20     50      100     150     200     250     280     300];
        M =     [0.46   0.64    0.78    0.68    0.44    0.23    0.18    0.18];
        uR2=abs(i2*R2);                                                  %abs(e-i1*R1-uC) // uR2=i2*R2


                                             %%% Select mtual inductance:
        Mn=0.8;
        %Mn = Poly_approximation(uL1, M, (uR2), 5);
        %Mn = Lagrange_Polynomial(uL1, M, (uR2));
        

        D1=L1/Mn-Mn/L2;
        D2=Mn/L1-L2/Mn;
        B = [-R1/(Mn*D1)   R2/(L2*D1)    -1/(Mn*D1);
             -R1/(L1*D2)  R2/(Mn*D2)     -1/(L1*D2);
                1/C             0           0   ];
        G = [1/(Mn*D1); 1/(L1*D2); 0];
        H = [i1; i2; uC];
        dy = B*H+G.*e;
    end


    function LaPoly = Lagrange_Polynomial(x, y, xp)
        LaPoly=0;
        z=x;     %vector of z's - nodes 
        kn = length(x);      %number of nodes
        for k=1:kn
            p=1;
            for j=1:kn
                if k~= j
                    p = p*(xp-z(j))/(z(k)-z(j));
                end
            end
            LaPoly =  LaPoly + y(k)*p;
        end
    end

    function p = Poly_approximation(x, y, xp, o)
           M=ones(length(x), o+1);
            for k=1:length(x)
                for j=2:o+1
                    M(k,j)=M(k,j-1)*x(k);
                end
            end
            %x = A\B solves the system of linear equations A*x = B. 
            Y = y';
            MTM = M \ M;
            MTY = M \ Y;
            A = MTM \ MTY;
                    
            A=(MTM)\(MTY);
            A=A';
            p=A(1);
            for k=2:o+1
                p=p+A(k)*xp^(k-1);
            end
        end
    
    function y1 = modified_euler(t, y1, h)
        for i=1:length(t)-1
            y1(:, i+1)=y1(:, i)+h*f(t(i)+h/2, y1(:, i)+f(t(i),y1(:, i))*h/2);
        end
    end

%integral - rectangle method

    function I = int_r(f, h) 
        l=length(f);
        I=zeros(1,l);
        for o=1:l
            I(1,o)=f(o)*h;
        end
    end

%integral - Simpson's method

    function I = int_p(f,h)
        l=length(f);
        I=zeros(1,l);
        for o=2:2:length(f)-1
            I(o)=2*(f(o))*h/3;
        end
        for o=3:2:length(f)-1
            I(o)=4*(f(o))*h/3;
        end
        I(1)=f(1)*h/3;
        I(end)=f(end)*h/3;
    end

    subplot(3,2,6);
        uR_2=zeros(4,length(t));
            uR_2(1,:)=R2.*i2_all(1,:);
            plot(t, uR_2(1,:));
            grid on;
            legend('UR2[V]');
            ylabel('UR2[V]');
            xlabel('t[s]');
            title('UR2(t)');

    for i=1:length(t)
        ev(i)=E(t(i));
    end
    subplot(3,2,1);
            plot(t,ev, '-');
            ylabel('e[V]');
            xlabel('t[s]');
            syms t;
            title(sprintf('e(t) = %s', E(t)));
end