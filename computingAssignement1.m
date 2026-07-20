%% 
% PART A
x = linspace(0.99, 1.01, 4001);

yf= ((x-1).^7);
ye = ((x.^7) - 7*(x.^6) + 21*(x.^5) - 35*(x.^4) + 35*(x.^3) -21*(x.^2) + 7*x - 1);
yh = ((((((x-7).*x + 21).*x -35).*x +35).*x -21).*x +7).*x -1;

figure("Name","Functions")
plot(x, ye, ".", x, yh, ".", x, yf, ".");
legend({"Pe(x)", "Ph(x)", "Pf(x)"});
title("Functions");
xlabel("X");
ylabel("Y");
clear;
% PART B
%% 
a = (linspace(0.95, 1.05, 1000));

yf= ((a-1).^7);
ye = ((a.^7) - 7*(a.^6) + 21*(a.^5) - 35*(a.^4) + 35*(a.^3) -21*(a.^2) + 7*a - 1);
yh = ((((((a-7).*a + 21).*a -35).*a +35).*a -21).*a +7).*a -1;

ee = abs(ye-yf);
eh = abs(yh-yf);

disp("Machine Precision: "+(eps))
figure("Name","Errors");
plot(abs(a-1), ee,".", abs(a-1), eh,".");
legend({"eE(x)", "eH(x)"});
title("Errors");
xlabel("|X-1|");
ylabel("Y");
clear


% Finding a valid delta
%% 
lastA = 0;
delta = -100;
for a = 0:0.0000001:2
    
    b = abs(a-1);
    Pe = abs((a^7) - 7*(a^6) + 21*(a^5) - 35*(a^4) + 35*(a^3) -21*(a^2) + (7*a) - 1);
    Pf = ((a-1)^7);
    res = abs(Pe - Pf)/abs(Pf);
    if ((res - 0.01 < (5*eps())))
        %a
        
        %disp("B: " + b);
        %disp("B-LastI: " + (b-lastI));
        %if(lastI ~= b-0.0001)
        if((abs((b-lastA)- 0.0000001) > (5*eps())))
            delta = b;
            %disp(newline+"LastI: " + lastI);
            %disp("B: " + b);
            %disp("B - 0.0001: " + (b-0.0000001));
        end
        lastA = b;
    end
end
disp("Delta: " + delta);
clear;


%For the delta section of Part B i need to design an matlab algorithms to
%find the right delta The delta is the distance from 1 which ensures that
%all values between 1 and 1+delta satisfy the error is greater than 0.01.
%So all of the errors outside are delta have error less than 0.01
