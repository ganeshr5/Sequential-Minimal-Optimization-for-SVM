
% SMO FOR SVM
% MNIST-13 data set

function [a] = mysmosvm(filename, numruns)

clear all
clc

filename='MNIST-13.csv';
numruns=10;

warning('OFF');

dat = csvread(filename);
n = size(dat,1);

ndat = dat(:, 2:end);  % new data except class labels
clabel = dat(:,1);     % class labels

for i = 1:n
    if (clabel(i) == 3)
        clabel(i) = 0;
    end
end


alpha2 = SMO.al(i2); y2 = SMO.y(i2);

if ((alpha2 > SMO.ep) & (alpha2 < (SMO.C-SMO.ep)))
    e2 = SMO.error(i2);
else
    e2 = fwd(i2) - y2;
end;

% r2 < 0 if point i2 is placed between margin (-1)-(+1)
% Otherwise r2 is > 0. r2 = f2*y2-1

r2 = e2*y2; 
%KKT conditions:

if ((r2 < -SMO.tolerance) & (alpha2 < (SMO.C-SMO.ep))) | ...
((r2 > SMO.tolerance) & (alpha2 > SMO.ep))
    % If it doens't violate KKT conditions then exit, otherwise continue.

    %Try i2 by three ways; if successful, then immediately return 1; 
    RESULT = 1;
    

    POS = find((SMO.al > SMO.ep) & (SMO.al < (SMO.C-SMO.ep)));
    [MAX,i1] = max(abs(e2 - SMO.error(POS)));
    if ~isempty(i1)
        if takeStep(i1, i2, e2), return;
        end;
    end;

    %The second heuristic choose any Lagrange Multiplier that is a SV and tries to optimize
    for i1 = randperm(SMO.ntp)
        if (SMO.al(i1) > SMO.ep) & (SMO.al(i1) < (SMO.C-SMO.ep))
        %if a good i1 is found, optimise
            if takeStep(i1, i2, e2), return;
            end;
        end
    end

    %if both heuristc above fail, iterate over all data set 
    for i1 = randperm(SMO.ntp)
        if ~((SMO.al(i1) > SMO.ep) & (SMO.al(i1) < (SMO.C-SMO.ep)))
            if takeStep(i1, i2, e2), return;
            end;
        end
    end;
end; 


a1 = alpha1 + s*(alpha2-a2);

w1 = y1*(a1 - alpha1); w2 = y2*(a2 - alpha2);



% To find objective function 
%objective function that depend on alpha2 need be evaluated... 

    ind = find(SMO.al>0);

    aa2 = L; aa1 = alpha1 + s*(alpha2-aa2);

    Lobj = aa1 + aa2 + sum((-y1*aa1/2).*SMO.y(ind).*K(ind,i1) + (-y2*aa2/2).*SMO.y(ind).*K(ind,i2));

    aa2 = H; aa1 = alpha1 + s*(alpha2-aa2);
    Hobj = aa1 + aa2 + sum((-y1*aa1/2).*SMO.y(ind).*K(ind,i1) + (-y2*aa2/2).*SMO.y(ind).*K(ind,i2));

    if (Lobj>Hobj+SMO.ep)
        a2 = H;
    elseif (Lobj<Hobj-SMO.ep)
        a2 = L;
    else
        a2 = alpha2;
    end;
end;

if (abs(a2-alpha2) < SMO.ep*(a2+alpha2+SMO.ep))
    return;
end;

a=5
end

end
