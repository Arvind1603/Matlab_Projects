%Name: Arvind Purohit




%-----------Extracted from D3-----------
clear


%Inputs parameters
R = 1;
while R == 1
    clear
    clf
    
    
    disp("Considering the differential equation: d^2f/dt^2 + g df/dt + hf(t) = acos(wt) + bsin(wt)")

    g = input('Input a value for g: '); 
    h = input('Input a value for h: ');  %Takes in the input
    a = input('Input a value for a: ');   
    b = input('Input a value for b: ');
    w = input('Input a value for w: ');
    
    %We form a matrix(2x2) 
    %we ended up with A(h-w^2) + Bgw = a and Agw + B(h-w^2) = b
    
    
    K_C = [a; b;]; % Known coefficents
    
    M = [h-w^2 g*w; -g*w h-w^2;]; %matrix format of the equation(inverse)
    
    U_C = M^(-1)*K_C;  %Unkwon coefficient

    %We established f(t) = Acos(wt)+Bsin(wt)
    A = U_C(1)   %coefficient of cos(wt)
    B = U_C(2)   %coefficient of sin(wt)
    
%---solved equation----

    tmax = 3*pi/w; 
    t = linspace(0,tmax,1000);
    
    f_t = A*cos(w*t) + B*sin(w*t); %Particular solution
    ff = a*cos(w*t) + b*sin(w*t);  %Forcing function
    
    amp_ff = max(ff);
    
    ratio = amp_ff/max(f_t);
    nscale = round(log10(ratio));
    
    ft = f_t*10^nscale;

    
    if tmax >= 1  % for the time parameter to be in second
        
        plot(t,ff,t,ft,'Linewidth',3)
        xlabel('Time t (s)','FontSize',18)
     
    elseif tmax < 1 && tmax > 0.001  % for the time parameter to be in millisecond
        tms = t*1e3;
        plot(tms,ff,tms,ft,'Linewidth',3)
        xlabel('Time t (ms)','FontSize',18)
        
    elseif tmax < 0.001 % for the time parameter to be in microsecond
        tus = t*1e6;
        plot(tus,ff,tus,ft,'Linewidth',3)
        xlabel('Time t (us)','FontSize',18)
    end
    
    if nscale == 0
        legend('forcing function','Response')
    else
        text = sprintf('Response (x10^%d)',nscale);
        legend('forcing function',text)
    end
    
    ax = gca; 
    ax.FontSize = 15;  %sets everything to be 15pt

    ylabel('f(t)','FontSize',18)
    
    ymax = 1.3*max(max(ff), max(ft));
    
    ylim([-ymax ymax])

    Text = sprintf('$ d^2f/dt^2 + %g df/dt + %gf(t) = %g \\cos(%gt) + %g \\sin(%gt) $',...
        g,h,a,w,b,w);
    
    title({'Plotting forcing function and response'},'FontSize',24)
    
    subtitle(Text,'Interpreter','latex','FontSize',24)
    
    R = input('Do you want to continue? (Press 1 for yes and 0 for no): ');

end
