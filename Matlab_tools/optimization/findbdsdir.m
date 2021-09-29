function [a,b,c,fail] = findbdsdir(fnct,a,b,dir,spot,minormax)

%  We define minormax as 1 or -1.  By using this, we can multiply any feval
%  by this constant (1 corresponds to a minimum, -1 to a maximum), in order
%  to find the minimum or maximum.

%  Assign some highly used numbers.
if a == b
    a = b+1;
end
gold = 2/(sqrt(5)-1);
glimit = 100;
zerobuffer = 10*eps;
bignum = 2^1023;
maxits = 1000;
count = 0;
exit = 0;
fail = 0;

%  Initialize function values of a and b.
fa = minormax*feval(fnct,a*dir+spot);
fb = minormax*feval(fnct,b*dir+spot);

%  Check to see if fa and fb are distinct, in other words is fa>fb or vice
%  versa?
if fa==fb
    
    %  If not, then search inside of interval for a point with a lower or
    %  higher value (search is not exhaustive but is extensive as in we try
    %  a lot of points via a sort of a binary slicing of the interval, so 
    %  it should find a point unless the fnct is weird, which, if the case,
    %  then the user is warned).
    fafbtest = 0;
    while fafbtest == 0
        middle = (b+a)/2;
        fm = minormax*feval(fnct,middle*dir+spot);
        if fm < fb
            c = b;
            fc = fb;
            b = middle;
            fb = fm;
            fafbtest = 1;
        elseif fm > fb
            b = a;
            fb = fa;
            a = middle;
            fa = fm;
            fafbtest = 1;
        elseif b-a > zerobuffer
            b = middle;
        else
%             error('Function seems to be constant by sampling of many points.');
            fail = 1;
            c = 0;
            return
        end
    end
    
    %  Make sure we're going downhill.  If not, switch and a and b and fa
    %  and fb.
elseif fb > fa
    temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
end

%  Start going downhill to find a c that takes us uphill.

c = b+gold*(b-a);
fc = minormax*feval(fnct,c*dir+spot);
while (fb >= fc) && (count < maxits)
    
    %  Find the u that corresponds to the bottom (hopefully) or top of the
    %  parabola given by a, b, and c.
    count = count + 1;
    r = (b-a)*(fb-fc);
    q = (b-c)*(fb-fa);
    if abs(q-r) > zerobuffer
        u = b-((b-c)*q-(b-a)*r)/(2*(q-r));
    else
        u = b-((b-c)*q-(b-a)*r)/(2*zerobuffer*sign(a-b));
    end
    fu = minormax*feval(fnct,u*dir+spot);
    
    %  Just define how far we'll go downhill before starting to think the
    %  whole thing is just one long downhill slope.
    ulim = b+glimit*(c-b);
    
    %  If u is between b and c, the bottom or top of the parab is in
    %  between the two.  If a bottom, store our new found min bounds.
    if ((b-u)*(u-c) > 0)
        if (fu < fc)
            a = b;
            fa = fb;
            b = u;
            fb = fu;
            exit = 1;
            
        %  If by some freak chance fu is actually above fb (even though if
        %  it really was a parabola asymptote, it would have to be the minimum) 
        %  then we might as well store this fact while we've already evaluated
        %  fu.
        elseif (fu < fb)
            c = u;
            fc = u;
            exit = 1;
        else

            %  Parabolic fit found no minimum bounds.  So, go downhill again:
            u = c+gold*(c-b);
            fu = minormax*feval(fnct,u*dir+spot);
        end
        
    %  u is beyond a, b, and c (but within proper bounds).  Check to see if
    %  it's not a minimum.  Shift b and c to the right it it's not.
    elseif ((c-u)*(u-ulim)>0)
        fu = minormax*feval(fnct,u*dir+spot);
        if fc > fu
            b = c;
            c = u;
            u = c+gold*(c-b);
            fb=fc;
            fc=fu;
            fu=minormax*feval(fnct,u*dir+spot);
        end
        
    %  We've gone past the edge!  What do we do?  Make the envelope bigger.
    elseif (u-ulim)*(ulim-c)>=0
        u = ulim;
        fu = minormax*feval(fnct,u*dir+spot);
        if glimit <= log(bignum)
            glimit = glimit*gold;
        else
            error('Can''t seem to find an increasing slope.');
            return;
        end
    
    %  If no previous results were used, then just step downhill again.
    else
        u = c+gold*(c-a);
        fu = minormax*feval(fnct,u*dir+spot);
    end
    
    %  Shift everything getting rid of the old value of a.
    if exit == 0
        a = b;
        b = c;
        c = u;
        fa = fb;
        fb = fc;
        fc = fu;
    end
    
    %  This is so we don't keep on stepping downhill "forever".
    if (abs(c)==Inf) || (abs(fc)==Inf)
        c = 'Too Big';
        return
    end
end
fb=minormax*fb;