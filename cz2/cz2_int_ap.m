function rlc
    clear;
    %close all;
    clc;

    F1=5;
    F2=50;

    syms t;
    E=@(t) 120*sin(2*pi*F2*t);
    %{ 
    E=@(t) 1;
    E=@(t) sin(t);
    E=@(t) 120*sin(2*pi*F2*t);
    E=@(t) 120*sin(2*t);
    E=@(t) 240*sin(2*t);
    %}

    h=0.001;
    a = 0;
    b = 30;
    t = a:h:b;      %[s]
    y=[0;           %[A]
       0;           %[A]   
       0];          %[V]
    R2 = 10;
    
                     %method number 1- poly interp, 2- spline interp; 3- 3rd order poly approx; 4- 5h order poly approx; 
                     %if m==5 -> M=const=0.8;
    for m=1:5
        %y=simply_euler(t, y, h);
        y1=modified_euler(t, y, h);                     %m                         
        for k=1:length(t)
            i1_all(m,:)=y1(1,:);
            i2_all(m,:)=y1(2,:);
            uC_all(m,:)=y1(3,:);
        end
    end

    function dy = f(t,y)
        e=E(t);
        R1 = 0.1;       %[Ohm]
        %R2 = 10;        %[Ohm]
        C = 0.5;        %[F]
        L1 = 3;         %[H]
        L2 = 5;         %[H]
        i1=y(1);
        i2=y(2);
        uC=y(3);
        uL1 =   [0      20     50      100     150     200     250     280     300];
        M =     [0.3    0.46   0.64    0.78    0.68    0.44    0.23    0.18    0.18];
        uR2=R2*i2;                                                 
        if     m==1
            Mn = Lagrange_Polynomial(uL1, M, (uR2));
        elseif m==2
            Mn = Spline_interp(uL1, M, (uR2));
        elseif m==3 
            Mn = Poly_approximation(uL1, M, (uR2), 3); 
        elseif m==4
            Mn = Poly_approximation(uL1, M, (uR2), 5);
        else
            Mn=interp1(uL1, M, uR2, "linear", "extrap");
        end
        D1=L1/Mn-Mn./L2;
        D2=Mn./L1-L2./Mn;
        B = [-R1/(Mn*D1)   R2/(L2*D1)    -1/(Mn*D1);
             -R1/(L1*D2)  R2/(Mn*D2)     -1/(L1*D2);
                1/C             0           0   ];
        G = [1/(Mn.*D1); 1/(L1.*D2); 0];
        H = [i1; i2; uC];
        dy = B*H+G.*e;
    end

    function LaPoly = Lagrange_Polynomial(x, y, xp)
        LaPoly=0;
        kn = length(x);
        for k=1:kn
            p=1;
            for j=1:kn
                if k~= j
                        if xp > x(end)
                            xp=x(end);
                        end
                        if xp < x(1)
                            xp=x(1);
                        end
                    p = p.*(xp-x(j))./(x(k)-x(j));
                end
            end
            LaPoly =  LaPoly + y(k)*p;
        end
end
    
    function S = Spline_interp(x, y, xp)        %1st order
        k=1;
        if xp > x(end-1)
            k=length(x)-1;
        else
            while xp >= x(k+1)
                k=k+1;
            end
        end
        A=[1 x(k);
           1 x(k+1)];
        B=[y(k);
           y(k+1)];
        X=A\B;
    
        S = X(1)+X(2)*xp;
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
        
        p=A(1);
        for k=2:o+1
            if xp > x(end)
                xp=x(end);
            end
            if xp < x(1)
                xp=x(1);
            end
            p=p+A(k).*xp.^(k-1);
        end
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
    
fig=figure;  
   fig.WindowState = 'maximized';
    subplot(3,2,3); 
        plot(t, i1_all(1,:), ...
            t, i1_all(2,:), ...
            t, i1_all(3,:), ...
            t, i1_all(4,:), ...
            t, i1_all(5,:));
        grid on;
        legend('Interpolacja metodą Lagrange`a', ...
            'Interpolacja splajnami 1 stopnia', ...
            'Aproksymacja wielomianami 3 stopnia', ...
            'Aproksymacja wielomianami 5 stopnia', ...
            'Interpolacja z ekstrapolacją, f. wbudowana, spline 1 st');
        ylabel('I[A]');
        xlabel('t[s]');
        title('I1(t)');
    
    subplot(3,2,4);
        plot(t, i2_all(1,:), ...
            t, i2_all(2,:), ...
            t, i2_all(3,:), ...
            t, i2_all(4,:), ...
            t, i2_all(5,:));
        grid on;
        legend('Interpolacja metodą Lagrange`a', ...
            'Interpolacja splajnami 1 stopnia', ...
            'Aproksymacja wielomianami 3 stopnia', ...
            'Aproksymacja wielomianami 5 stopnia', ...
            'Interpolacja z ekstrapolacją, f. wbudowana, spline 1 st');
        ylabel('I[A]');
        xlabel('t[s]');
        title('I2(t)');
        
    subplot(3,2,5);
        plot(t, uC_all(1,:), ...
            t, uC_all(2,:), ...
            t, uC_all(3,:), ...
            t, uC_all(4,:), ...
            t, uC_all(5,:));
        grid on;
        legend('Interpolacja metodą Lagrange`a', ...
            'Interpolacja splajnami 1 stopnia', ...
            'Aproksymacja wielomianami 3 stopnia', ...
            'Aproksymacja wielomianami 5 stopnia', ...
            'Interpolacja z ekstrapolacją, f. wbudowana, spline 1 st');
        ylabel('UC[V]');
        xlabel('t[s]');
        title('UC(t)');

    subplot(3,2,6);
            uR_2=zeros(5,length(t));
            for g=1:5
                uR_2(g,:)=R2.*i2_all(g,:);
            end
            plot(t, uR_2(1,:), ...
                t, uR_2(2,:), ...
                t, uR_2(3,:), ...
                t, uR_2(4,:), ...
                t, uR_2(5,:));
            grid on;
            legend('Interpolacja metodą Lagrange`a', ...
                'Interpolacja splajnami 1 stopnia', ...
                'Aproksymacja wielomianami 3 stopnia', ...
                'Aproksymacja wielomianami 5 stopnia', ...
                'Interpolacja z ekstrapolacją, f. wbudowana, spline 1 st');
            ylabel('UR2[V]');
            xlabel('t[s]');
            title('UR2(t)');

    subplot(3,2,1);
            plot(t,ev, '-');
            ylabel('e[V]');
            xlabel('t[s]');
            syms t;
            title(sprintf('e(t) = %s', E(t)));

    subplot(3,2,2);
        ylabel('M');
        xlabel('uL');
        title('M(uL)')
            uL1k =   [0      20     50      100     150     200     250     280     300];
            Mk =     [0.3    0.46   0.64    0.78    0.68    0.44    0.23    0.18    0.18];
            uL= 0 : 1 : uL1k(end);
                for i=1:length(uL);
                    Mia(1,i) = Lagrange_Polynomial(uL1k, Mk, uL(i));
                    Mia(2,i) = Spline_interp(uL1k, Mk, uL(i));
                    Mia(3,i) = Poly_approximation(uL1k, Mk, uL(i), 3);
                    Mia(4,i) = Poly_approximation(uL1k, Mk, uL(i), 5);
                    Mia(5,i) = interp1(uL1k, Mk, uL(i), 'linear');
                end
            hold on;
            plot(uL, Mia(1,:), ...
                uL, Mia(2,:), ...
                uL, Mia(3,:), ...
                uL, Mia(4,:), ... 
                uL, Mia(5,:), '--', ...
                uL1k, Mk, 'o');
                grid on;
            legend('Interpolacja metodą Lagrange`a', ...
                'Interpolacja splajnami 1 stopnia', ...
                'Aproksymacja wielomianami 3 stopnia', ...
                'Aproksymacja wielomianami 5 stopnia', ...
                'Interpolacja z ekstrapolacją, f. wbudowana, spline 1 st');
end