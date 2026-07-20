format longG;
function bijectionFunc = bisectionFunc(x)
    arguments
        x double;
    end
    bijectionFunc = x - cos(x);
    return;
end
function derivativeFunc = derivativeFunc(x)
    arguments
        x double;
    end
    derivativeFunc = 1 + sin(x);
    return;
end
function [estimate, iterations] = bisectionPrep(func, interval)  
    %disp(interval(1));
    %disp(interval(2));
    %disp("Prep func: " +func(1));
    if((interval(1) > interval(2)))
        warning("Invalid Interval");
        estimate = inteval(1) +1 ;
        iterations =[];
        return;
    end
    if((func(interval(1)) * func(interval(2))) > 0)
        warning("No Roots");
        estimate = interval(2)+1;
        iterations =[];
        return;
    end
    [estimate, iterations] = bisection(func, interval, 1);
    return;
end
function [estimate, iterations] = bisection(func, interval, currentIteration)
   mid = 0.5*(interval(1)+interval(2));
   if(currentIteration == 100 || (0.5*(interval(2)-interval(1))) < (10^-12))
      estimate = mid;
      iterations(currentIteration) = mid;
      disp("Bisection Total Iterations: " + (currentIteration));
      return;
   end
   fmid = func(mid);
   fint1 = func(interval(1));
   if((fint1*fmid) < 0)
       [estimate, iterations] = bisection(func, [interval(1), mid], currentIteration+1);
       iterations(currentIteration) = mid;
       return;
   else
       [estimate, iterations] = bisection(func, [mid, interval(2)], currentIteration+1);
       iterations(currentIteration) = mid;
       return;
   end
end
function [estimate, iterations] = fixedPoint(func, curX)
    prevX = -1;
    currentIteration = 1;
    while (abs(curX - prevX) >= (10^-12) && currentIteration < 100)
        iterations(currentIteration) = curX; 
        prevX = curX;
        curX = func(curX);
        currentIteration = currentIteration+1;
    end
    disp("Fixed Point Iterations: "+ (currentIteration-1));
    iterations(currentIteration) = curX;
    estimate = curX;
    return;
end
function [estimate, iterations] = newtonsMethod(func, deriv, curX)
    currentIteration = 1;
    prevX = 1;
    while(abs(curX-prevX) >= (10^-12) && currentIteration < 100)
        iterations(currentIteration) = curX;
        prevX = curX;
        curX = curX - (func(curX)/deriv(curX));
        currentIteration = currentIteration +1;
    end
    iterations(currentIteration) = curX;
    estimate = curX;
    disp("Newtons Method Total Iterations: " + (currentIteration-1));
    return;
end
function [estimate, iterations] = secantMethod(func, curX, prevX)
    currentIteration = 1;
    while(abs(curX-prevX) >= (10^-12) && currentIteration < 100)
        iterations(currentIteration) = curX;
        newX = curX - (func(curX)*((curX-prevX)/(func(curX)-func(prevX))));
        prevX = curX;
        curX = newX;
        currentIteration = currentIteration +1;
    end
    iterations(currentIteration) = curX;
    estimate = curX;
    disp("Secant Method Total Iterations: " + (currentIteration - 1));
    return;
end
[bEstimate,bIterations] = bisectionPrep(@bisectionFunc, [0, 1]);
[fEstimate, fIterations] = fixedPoint(@cos, 0.5);
[nEstimate, nIterations] = newtonsMethod(@bisectionFunc,@derivativeFunc, 0.5);
[sEstimate, sIterations] = secantMethod(@bisectionFunc, 1, 0);
%Both newtons and secant narrow down on x* to the same level of accuracy so
%im going to be using newtons method but secant would work too.
mostAccurate = nEstimate
figure("Name","Erros For Estimates of Cos(x)");
semilogy(abs(bIterations-mostAccurate), "b");
hold on;
semilogy(abs(fIterations-mostAccurate), "r");
semilogy(abs(nIterations-mostAccurate), "g");
semilogy(abs(sIterations-mostAccurate), "y");
title("Errors");
legend({"Bisection", "Fixed Point", "Newtons", "Secant"});
xlabel("N");
ylabel("Accuracy");
hold off;
grid("on");
function [x,y] = getArraysForPolyfit(iterations, mostAccurate, offset)
    initial = floor(0.5*length(iterations))-1;
    y = zeros(length(iterations)-offset-initial);
    x = zeros(length(iterations)-offset-initial);
    for i = initial:length(iterations)-offset-1
        y(i-initial+1) = log(abs(iterations(i+1)-mostAccurate));
        x(i-initial+1) = log(abs(iterations(i)-mostAccurate));
    end
end
[fX,fY] = getArraysForPolyfit(fIterations, mostAccurate, 0);
fEstimate = polyfit(fX, fY, 1);
disp("Fixed Point");
disp("Alpha Estimate");
fEstimate(1)
disp("Lambda Estimate");
exp(fEstimate(2))
[bX,bY] = getArraysForPolyfit(bIterations, mostAccurate, 0);
bEstimate = polyfit(bX, bY, 1);
disp("Bisection");
disp("Alpha Estimate");
bEstimate(1)
disp("Lambda Estimate");
exp(bEstimate(2))
[nX,nY] = getArraysForPolyfit(nIterations, mostAccurate, 2);
nEstimate = polyfit(nX, nY, 1);
disp("Newtons");
disp("Alpha Estimate");
nEstimate(1)
disp("Lambda Estimate");
exp(nEstimate(2))
[sX,sY] = getArraysForPolyfit(sIterations, mostAccurate, 1);
sEstimate = polyfit(sX, sY, 1);
disp("Alpha Estimate");
sEstimate(1)
disp("Lambda Estimate");
exp(sEstimate(2))
clear;


% OLD TESTING FUNCTIONS IF I EVER NEED THEM AGAIN.

function lambda = findLambda(Ekplusone, Ek,mostAccurate, alpha)
lambda = abs(Ekplusone-mostAccurate)/(abs(Ek-mostAccurate)^alpha);
end

function [x,y] = testGetArraysForPolyfit(iterations, mostAccurate, offset)
    for j = 1:1:(length(iterations)-offset)
        %initial = floor(0.5*length(iterations))-1;
        initial = j;
        y = zeros(length(iterations)-offset-initial);
        x = zeros(length(iterations)-offset-initial);

        for i = initial:1:length(iterations)-offset-1
            y(i-initial+1) = log(abs(iterations(i+1)-mostAccurate));
            x(i-initial+1) = log(abs(iterations(i)-mostAccurate));
        end
        estimate = polyfit(x, y, 1);
        disp("Alpha: " + estimate(1));
        disp("Lambda: " + exp(estimate(2)));
    end
end

%Testing functions so i can understand what happends to my estimatesion as
%I move to different starting points. From what i see my data consistently
%gets a simmilar value with starting form the very begining being slightly worse
%so im going to stick with my simple floor(length/2) -1 starting point to
%avoid having my estimate thrown off by the less accurate low values.
%disp("TESTING DIFFERENT STARTING POINTS FOR FIXED POINT");
%testGetArraysForPolyfit(fIterations, mostAccurate, 0);
%disp("TESTING DIFFERENT STARTING POINTS FOR BISECTION");
%testGetArraysForPolyfit(bIterations, mostAccurate, 0);
%disp("TESTING DIFFERENT STARTING POINTS FOR NEWTONS");
%testGetArraysForPolyfit(nIterations, mostAccurate, 2);
%disp("TESTING DIFFERENT STARTING POINTS FOR NEWTONS");
%testGetArraysForPolyfit(sIterations, mostAccurate, 1);


% I NEED TO FIGURE OUT WHY MY LAMBDAS ARE SO BAD. MAYBE GO BACK TO THE OLD
% METHOD OF TESTING SOME HIGH K LAMBDAS AND SEEING WHERE THEY GENREALLY
% CONVERGED FOR VALUES NEAR MY ESTIMATE

function testForLambda(iterations, offset, alphaEstimate, mostAccurate)
for i = alphaEstimate-0.01:0.00001:alphaEstimate+0.01
    for j = 0:1:2
        lambda = abs(iterations((end-offset)-j)-mostAccurate)/(abs(iterations(((end-offset)-j)-1)-mostAccurate)^i);
        disp("Lambda: "+ lambda);
    end
end
end

%clear;


