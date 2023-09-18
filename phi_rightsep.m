function y = phi_rightsep(i,j,n1,n2,nf,x)
N2 = n1+n2;
h1 = nf/(n1*(n1+nf));
low = max(0,h1-(j/N2-i/n1));
knot_original = [(i-1)/n1,(i-1)/n1+h1-low,j/N2-1/(n1+nf)-(h1-low),j/N2-1/(n1+nf),(j-1)/N2];
bottom_original = [h1,low,low,h1,h1];
for oo = 2:5
    if (j-2)/N2 >= knot_original(oo-1) && (j-2)/N2<knot_original(oo)
        knot = [(j-2)/N2,knot_original(oo:5)];
        b = bottom_original(oo-1)+((j-2)/N2-knot_original(oo-1))*(bottom_original(oo)-bottom_original(oo-1))/(knot_original(oo)-knot_original(oo-1));
        bottom = [b,bottom_original(oo:5)];
    end
end
t = length(knot);
y = 0;
for oo = 1:t-1
    if x > knot(oo+1)
        y = y+(knot(oo+1)-knot(oo))/2*(bottom(oo+1)+bottom(oo));
    elseif knot(oo+1)-knot(oo) ~= 0
        y = y+(x-knot(oo))/2*(bottom(oo)+(bottom(oo+1)-bottom(oo))/(knot(oo+1)-knot(oo))*(x-knot(oo))+bottom(oo));
        break;
    end
end
end