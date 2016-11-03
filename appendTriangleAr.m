function T = appendTriangleAr(T, trAr)
global boundCond
global xL xR yL yR
global trCnt

for i=1:8
    if (trAr(i,2) < xL) || (trAr(i,4) < xL) || (trAr(i,6) < xL) 
        boundCond(1) = 1;
    elseif (trAr(i,2) > xR) || (trAr(i,4) > xR) || (trAr(i,6) > xR)
        boundCond(2) = 1;
    elseif (trAr(i,3) < yL) || (trAr(i,5) < yL) || (trAr(i,7) < yL)
        boundCond(3) = 1;
    elseif (trAr(i,3) > yR) || (trAr(i,5) > yR) || (trAr(i,7) > yR)
        boundCond(4) = 1;
    else
        T = [T; trAr(i,:)];
        trCnt = trCnt + 1;
    end
    
end