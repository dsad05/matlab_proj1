function rlc
    clear;
    %close all;
    clc;
tic
    aF=0.0;
    bF=1;
    hF=0.1;
    F=aF:hF:bF;
    p_0=406;

    for io=1:length(F)
        syms t;
        E=@(t) 100*sin(2*pi*F(io)*t);

        epsx=1e-6;
        epsy=1e-6;

        h=0.001;
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
        P_inst=zeros(1,length(t));
        y=modified_euler(t, y, h, E);                      
        for k=1:length(t)
            P_inst(1,:)=(y(1,:).^2).*R1+(y(2,:).^2).*R2;
        end
    
        I_p = int_p(P_inst, h);
        I_p_sum=sum(I_p);
        Pp(1,io)=I_p_sum;
    end

    P=Pp-p_0;
    interv=F;
    df=0.001;
    toc
%BISECTION
    tic
    fprintf('\tBisection method:\n')
    x_0 = bisection(interv, P, epsx, epsy, y, h, p_0);

    fprintf('\tFor interval <%.1f ; %.1f> ', aF,bF)
    if (isempty(x_0))
        fprintf('there is no roots\n')
    elseif (length(x_0)==1)
        fprintf('there is %.f root:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    else
        fprintf('there are %.f roots:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    end
 fprintf('___________________________________________________________________\n\n')
    toc

%SECANT
    tic
    fprintf('\tSecant method:\n')
    x_0 = secant(interv, P, epsx, epsy, y, h, p_0);
    fprintf('\tFor interval <%.1f ; %.1f> ', aF,bF)
    if (isempty(x_0))
        fprintf('there is no roots\n')
    elseif (length(x_0)==1)
        fprintf('there is %.f root:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    else
        fprintf('there are %.f roots:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    end
     fprintf('___________________________________________________________________\n\n')
    toc

%QUASI-NEWTON
    tic
    fprintf('\tQuasi-Newton method:\n')
    x_0 = q_newton(interv, P, epsx, epsy,  y, h, p_0);
    fprintf('\tFor interval <%.1f ; %.1f> ', aF,bF)
    if (isempty(x_0))
        fprintf('there is no roots\n')
    elseif (length(x_0)==1)
        fprintf('there is %.f root:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    else
        fprintf('there are %.f roots:\n', length(x_0))
        fprintf('\tx_0 = %.4f\n', x_0);
    end
     fprintf('___________________________________________________________________\n\n')
    toc
        fig=figure;  
    fig.WindowState = 'maximized';
        subplot(1,1,1);         
            plot(F, Pp);
            yline(p_0);
            grid on;
            legend('P');
            ylabel('P[W]');
            xlabel('f[Hz]');
            title('P(f)');

%FUNCTIONS
    function f_diff = fdiff(f, df,h,y,p_0)
        f_diff = (P0(f+df,h,y,p_0)-P0(f,h,y,p_0))/df;
    end

    function ipsum = P0(x_i, h, y, p_0)
        syms t;
        E=@(t) 100*sin(2*pi*x_i*t); 
        a = 0;
        b = 30;
        t = a:h:b;      %[s]
        E(t);
        y1=modified_euler(t, y, h, E);                      
        Pinst(1,:)=(y1(1,:).^2).*R1+(y1(2,:).^2).*R2;
        Ip = int_p(Pinst, h);
        ipsum=sum(Ip)-p_0;
    end

    function x_0 = bisection(interv, p, epsx, epsy, y, h, p_0)
        x_0=[];
        Pcounter=0;
        bcounter=0;
            for i=1:length(interv)-1
                ai=interv(i);
                bi=interv(i+1);
                fa=p(i);
                fb=p(i+1);
                fprintf('Interval:')
                fprintf('\t %.2f', ai, bi)
                fprintf('\n')
                if (fa * fb > 0)
                    disp('The same sign at endpoints of the interval')
                else
                    counter=0;
                    while(true)
                        counter=counter+1;
                        x_i=(ai+bi)/2;
                        if abs(ai-x_i) < epsx
                            break;
                        end
    
                        f_div=P0(x_i, h, y, p_0);
                        Pcounter=Pcounter+1;
    
                        if abs(f_div) < epsy
                            break;
                        end
                        if (fa*f_div < 0)
                            bi = x_i;
                        else
                            ai = x_i;
                            fa=f_div;
                        end
                    end
                    bcounter=bcounter+counter;
                    x_0(end+1)=x_i;
                    fprintf('Root:\t\t')
                    fprintf('%.5f \n',x_i)
                    fprintf('Iteration counter per root: %.f \n',counter)
                end
                fprintf('________________________________________________\n\n')
            end
        fprintf('P-function counter per bisection function: %.f \n',Pcounter)
        fprintf('Iteration counter per bisection function: %.f \n',bcounter)
    end
    
    function x_0 = secant(interv, p, epsx, epsy, y, h, p_0)
        x_0=[];
        Pcounter=0;
        bcounter=0;
        for i=1:length(interv)-1
            an=interv(i);
            bn=interv(i+1);
            fa=p(i);
            fb=p(i+1);
            fprintf('Interval:')
            fprintf('\t %.2f', an, bn)
            fprintf('\n')
            if (fa * fb > 0)
                disp('The same sign at endpoints of the interval')
            else
                if fa*(fdiff(fdiff(an, df,h,y,p_0),df,h,y,p_0)) > 0
                    Pcounter=Pcounter+4;
                    x_i = an;
                    fi=fa;
                else
                    x_i = bn;
                    fi=fb;
                end
                counter=0;
                while(true)
                    counter=counter+1;
                    if x_i == 0
                        x_i = 0.01;
                    end
                    x_h=x_i-df;
                    x_j=x_i-(fi*(x_i-x_h)/(P0(x_i, h, y, p_0)-P0(x_h, h, y, p_0)));
                    fj=P0(x_j, h, y, p_0);
                    Pcounter=Pcounter+3;
                    if abs(x_j-x_i) < epsx
                        break;
                    end
                    if abs(fj) < epsy
                        break;
                    else
                        x_i=x_j;
                        fi=P0(x_i, h, y, p_0);
                        Pcounter=Pcounter+1;
                    end
                end
                x_0(end+1)=x_i;
                fprintf('Root:\t\t')
                fprintf('%.5f \n',x_i)
                fprintf('Iteration counter per root: %.f \n',counter)
                bcounter=bcounter+counter;
            end
            fprintf('_____________________________________________________\n\n')
        end
        fprintf('P-function counter per secant function: %.f \n',Pcounter)
        fprintf('Iteration counter per secant function: %.f \n',bcounter)
    end

    function x_0 = q_newton(interv, p, epsx, epsy,  y, h, p_0)
    x_0=[];
    Pcounter=0;
    bcounter=0;
    counter = 0;
        for i=1:length(interv)-1
            an=interv(i);
            bn=interv(i+1);
            fa=p(i);
            fb=p(i+1);
            fprintf('Interval:')
            fprintf('\t %.2f', an, bn)
            fprintf('\n')
            if (fa * fb > 0)
                disp('The same sign at endpoints of the interval')
            else
                if fa*(fdiff(fdiff(an, df,h,y,p_0),df,h,y,p_0)) > 0
                    Pcounter=Pcounter+4;
                    x_i = an;
                    fi=fa;
                else
                    x_i = bn;
                    fi=fb;
                end
                counter=0;
                while(true)
                    counter=counter+1;
                    if x_i == 0
                        x_i = 0.01;
                    end
                    x_j=x_i-(fi/fdiff(x_i, df,h,y,p_0));
                    Pcounter=Pcounter+2;
                    fj=P0(x_j, h, y, p_0);
                    Pcounter=Pcounter+1;
                    if abs(x_j-x_i) < epsx
                        break;
                    end
                    if abs(fj) < epsy
                        break;
                    else
                        x_i=x_j;
                        fi=P0(x_i, h, y, p_0);
                        Pcounter=Pcounter+1;
                    end
                end
                bcounter=bcounter+counter;
                x_0(end+1)=x_i;
                fprintf('Root:\t\t')
                fprintf('%.5f \n',x_i)
                fprintf('Iteration counter per root: %.f \n',counter)
            end
            fprintf('_____________________________________________________\n\n')
        end
        fprintf('P-function counter per quasi-Newton function: %.f \n',Pcounter)
        fprintf('Iteration counter per quasi-Newton function: %.f \n',bcounter)
    end
    
    function dy = f(t,y,e)
        i1=y(1);
        i2=y(2);
        uC=y(3);
        uR2=R2*i2;                                               
        Mn=0.8;
        D1=L1/Mn-Mn./L2;
        D2=Mn./L1-L2./Mn;
        B = [-R1/(Mn*D1)   R2/(L2*D1)    -1/(Mn*D1);
             -R1/(L1*D2)  R2/(Mn*D2)     -1/(L1*D2);
                1/C             0           0   ];
        G = [1/(Mn.*D1); 1/(L1.*D2); 0];
        H = [i1; i2; uC];
        dy = B*H+G.*e;
    end

    function y1 = modified_euler(t, y1, h, E)
        for i=1:length(t)-1
            y1(:, i+1)=y1(:, i)+h*f(t(i)+h/2, y1(:, i), E(t(i)+h/2)+f(t(i),y1(:, i), E(t(i)))*h/2);
        end
    end

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

end
