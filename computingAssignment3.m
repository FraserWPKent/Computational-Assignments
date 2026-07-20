%clc
%% 
format long;
denseN = [100, 200, 400, 800, 1600, 3200];
timingArray = [1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3;];
for n = 1:1:length(denseN)
    for i = 1:1:5
        R = randn(denseN(n));
        A = R + denseN(n)*eye(denseN(n));
        b = rand(denseN(n),1);
        tic;
        inv(A)*b;
        curToc = toc;
        if(curToc < timingArray(n,1))
            timingArray(n,1) = curToc;
        end
        tic;
        A\b;
        curToc = toc;
        if(curToc < timingArray(n,2))
            timingArray(n,2) = curToc;
        end
    end
end
sparseN = [1000, 3000, 10000, 30000, 100000, 300000];
for n = 1:1:length(sparseN)
    for i = 1:1:5
        e = ones(sparseN(n),1);
        A = spdiags([e, -2*e, e], -1:1, sparseN(n), sparseN(n));
        b=randn(sparseN(n),1);
        tic;
        A\b;
        curToc = toc;
        if(curToc < timingArray(n,3))
            timingArray(n,3) = curToc;
        end
    end
end
method1 = zeros(1);
method2 = zeros(1);
method3 = zeros(1);

for i = 1:1:6
   method1(i)=timingArray(i,1);
   method2(i)=timingArray(i,2);
   method3(i)=timingArray(i,3);
end
figure("Name","Part A Timings");
loglog(method1, "r");
hold on;
loglog(method2, "g");
loglog(method3, "b");
hold off;
xlabel("N");
ylabel("Time");
xticks([1,2,3,4,5,6]);
xtickangle(15);
ticks = ["100/1000", "200/3000", "400/10000", "800/30000", "1600/100000","3200/300000"];
xticklabels(ticks);
title("Timings");
legend("Inv","Dense Backslash","Sparse Backslash","Location","northwest");

invEstimate = (log10(method1(6))-log10(method1(1)))/(log10(1600)-log10(100))
denseBackslashEstimate = (log10(method2(6))-log10(method2(1)))/(log10(1600)-log10(100))
sparseBackslashEstimate = (log10(method3(6))-log10(method3(1)))/(log10(300000)-log10(1000))

%%
mValues = [5, 20, 50, 100, 200];
timingArrayB = [5 5; 5 5; 5 5; 5 5; 5 5;];
n = 800;
for x = 1:1:length(mValues)
    for i = 1:1:5
        R = randn(n);
        A = R + n*eye(n);
        totalTime = 0;
        for j = 1:1:mValues(x)
            b_j = rand(n,1);
            tic;
            A\b_j;
            totalTime = totalTime + toc;
        end
        if(i==1 || totalTime < timingArrayB(x,1))
            timingArrayB(x,1) = totalTime;
        end
        totalTime = 0;
        [L,U,P] = lu(A);
        for j = 1:1:mValues(x)
            b_j = rand(n,1);
            tic;
            U\(L\(P*b_j));
            totalTime = totalTime + toc;
        end
        if(i==1||totalTime < timingArrayB(x,2))
            timingArrayB(x,2) = totalTime;
        end 
    end
end
method1B = zeros(1);
method2B = zeros(1);
for i = 1:1:5
    method1B(i)=timingArrayB(i,1);
    method2B(i)=timingArrayB(i,2);
end
figure("Name","Part B Timings");
loglog(method1B, "r");
hold on;
loglog(method2B, "g");
hold off;
title("Timings");
xlabel("M");
ylabel("Time");
xticks([1,2,3,4,5]);
xticklabels({"5", "20", "50", "100", "500"});
legend("Naive","Reuse","Location","northwest");
disp("Speed Difference Between Naive and Reuse:");
disp("Naive Time: ")
method1B(5)
disp("Reuse Time: ")
method2B(5)
disp("Difference Between the two")
method1B(5)-method2B(5)
clear

%slopeMethod1 = [method1(5), method1(6)];
%i=[5,6];
%polyfit(slopeMethod1,i, 1)
%slopeMethod2 = [method2(5), method2(6)];
%polyfit(slopeMethod2,i, 1)
%slopeMethod3 = [method3(5), method3(6)];
%%polyfit(slopeMethod3,i, 1)

%slope1 = 10^(((method1(6)-method1(1)))/(6-1))
%slope2 = 10^(((method2(6)-method2(1)))/(6-1))
%slope3 = 10^(((method3(6)-method3(1)))/(6-1))
%slope2 = log((method1(6)-method1(1))/(denseN(6)-denseN(1)))

%Holder for old wording from word doc
%As every time our value m increases by k values the naive implementation will increase by k(n^3) operations whereas the reuse implementation will only increase by k(n^2). This means that assuming n is constant the rate that our naive implementation grows will be n times large than the reuse implementations. This means that every time the value m goes up the difference between naive and reuse will grow significantly which in turn leads to a growing speed difference between the two


