function rlc
    clear;
    %close all;
    clc;

    T=3;
    F1=5;
    F2=50;

    syms t;
    E=@(t) 1;
    %{ 
    E=@(t) 1;
    E=@(t) sin(t);
    E=@(t) period(t,T);
    E=@(t) 240*sin(t);
    E=@(t) 210*sin(2*pi*F1*t);
    E=@(t) 120*sin(2*pi*F2*t);
    %}
    h=0.001;
    a = 0;
    b = 30;
    t = a:h:b;      %[s]
    y=[0;           %[A]
       0;           %[A]   
       0];          %[V]

    function e = period(t,T)
    p=rem(t,T);
        if(p>= T/2)
            e=0;
        else
            e=120;
        end
    end

    y=simply_euler(t, y, h);
    y1=modified_euler(t, y, h);

   fig=figure;  
   fig.WindowState = 'maximized';
    subplot(2,2,3);                     %simply euler i1 i2
        plot(t, y(1,:), '--', ...
            t, y(2,:), ...
            t, y(3,:))          %simply euler uC
        grid on;
        legend('I1[A]','I2[A]', 'UC[V]');
        ylabel('I[A], U[V]');
        xlabel('t[s]');
        title('I1(t), I2(t), UC(t) - prosta metoda Eulera');
    
    subplot(2,2,4);                      
        plot(t, y1(1,:), '--', ...
            t, y1(2,:), ...
            t, y1(3,:)) ;        %modified euler uC
        grid on;
        legend('I1[A]','I2[A]', 'UC[V]');
        ylabel('I[A], U[V]');
        xlabel('t[s]');
        title('I1(t), I2(t), UC(t) - zmodyfikowana metoda Eulera');
    function dy = f(t,y)
        e=E(t);
        R1 = 0.1;       %[Ohm]
        R2 = 10;        %[Ohm]
        C = 0.5;        %[F]
        L1 = 3;         %[H]
        L2 = 5;         %[H]
        M = 0.8;        %[H]
        i1=y(1);
        i2=y(2);
        uC=y(3);
        D1=L1/M-M/L2;
        D2=M/L1-L2/M;
        
        B = [-R1/(M*D1)   R2/(L2*D1)    -1/(M*D1);
             -R1/(L1*D2)  R2/(M*D2)     -1/(L1*D2);
                1/C             0           0   ];
        G = [1/(M*D1); 1/(L1*D2); 0];
        H = [i1; i2; uC];
        dy = B*H+G.*e;
    end

    function y = simply_euler(t, y, h)
        for i=1:length(t)-1
            y(:, i+1)=y(:, i)+h*f(t(i),y(:, i));
        end
    end
    
    function y1 = modified_euler(t, y1, h)
        for i=1:length(t)-1
            y1(:, i+1)=y1(:, i)+h*f(t(i)+h/2, y1(:, i)+f(t(i),y1(:, i))*h/2);
        end
    end

    for i=1:length(t)
        ev(i)=E(t(i));
    end
    subplot(2,2,1.5);    
        plot(t,ev, '-');
        ylabel('e[V]');
        xlabel('t[s]');
        title('e(t)');
        %ylim([-5 125]);
        syms t;
        title(sprintf('e(t) = %s', E(t)));
end