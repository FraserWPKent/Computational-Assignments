format("longG");
getJacobisWithC(2, "Part A 200");
clear;
%%
getJacobisWithC(4, "Part B 200");
function getJacobisWithC(c, name)
    iterationsArray = zeros(0,1);
    startNValues = 10;
    endNValues = 200;
    for n = startNValues:10:endNValues
        h = 1 / (n + 1);
        main_diag = (c/h^2) * ones(n,1);
        off_diag = (-1/h^2) * ones(n,1);
        A_n = spdiags([off_diag main_diag off_diag], [-1 0 1], n, n);
        b_n = ones(n, 1);
        xExact = (A_n)\b_n;
        iterations = (jacobiMethod(A_n, main_diag, b_n, n, xExact));
        iterationsArray((n/10)) = iterations;
    end
    figure("Name",name);
    loglog(iterationsArray, "r");
    if(c==2)
        hold on;
        xSquared = zeros(0,1);
        xCubed = zeros(0, 1);
        for n = startNValues:10:endNValues
            xSquared(end+1) = n^2;
            xCubed(end+1) = n^3;
        end
        loglog(xSquared, "g");
        loglog(xCubed, "b");
        hold off;
        legend("Number Of Iterations","N^2","N^3");
    else
        legend("Number Of Iterations");
    end
    legend("Location", "northwest")
    xlabel("N");
    ylabel("Iterations");
    title("C = " + c);
    grid;
    xticks(1:1:length(iterationsArray));
    xtickangle(90);
    tickValues = startNValues:10:endNValues;
    for i = 1:1:length(iterationsArray)
        tickValues(i) = string(tickValues(i));
    end
    xticklabels(tickValues);
    (polyfit((tickValues), iterationsArray, 1))
end

function iterations = jacobiMethod(matrixA, main_diag, b, n, xExact)
    matrixAWithoutD = matrixA-diag(main_diag);
    inverse = inv(diag(main_diag));
    iterations = 0;
    xK = zeros(n,1);
    while true
        xK = inverse*(b-(matrixAWithoutD * xK));
        xNum = norm((xExact - xK), inf);
        xDenom = norm(xExact,inf);
        xDiff = xNum/xDenom;
        iterations = iterations + 1;
        if xDiff < (10^-3) || iterations >= 1000000
            break;
        end
    end
end
clear;