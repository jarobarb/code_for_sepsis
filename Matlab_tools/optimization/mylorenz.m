function myf = mylorenz(x,Ab,At,x0,w)
    myf = Ab+(At-Ab)/(1+exp(-(x-x0)/w));