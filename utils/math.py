import decimal as dp

pi = dp.Decimal('3.141592653589793238462643383279')

###
###    trigonometric functions
###

def sin_d(x):
    """
    Returns the sine of x as measured in radians.
    """

    dp.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1

    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign

    dp.getcontext().prec -= 2

    return +s

def arcsin_d(x):
    """
    Returns the inverse sine of x as measured in radians.
    """

    dp.getcontext().prec += 2
    i, lasts, s, fact, num = 1, 0, x, 1, x

    while (s != lasts):
        lasts = s
        i += 2
        fact *= dp.Decimal(((i - 1) / 2))**dp.Decimal('2') / ((i - 1) * (i - 2))
        num *= x * x
        s += num / fact / i / 4**dp.Decimal((i - 1) / 2)

    dp.getcontext().prec -= 2

    return +s

def cos_d(x):
    """
    Returns the cosine of x as measured in radians.
    """

    dp.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1

    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign

    dp.getcontext().prec -= 2

    return +s

def arccos_d(x):
    """
    Returns the inverse sine of x as measured in radians.
    """
    
    dp.getcontext().prec += 2
    ac = pi / dp.Decimal('2') - arcsin_d(x)
    dp.getcontext().prec -= 2

    return ac

def tan_d(x): 
    """
    Returns the tangent of x as measured in radians.
    """

    dp.getcontext().prec += 2
    i, lasts, s, fact, num = 0, 0, x, 1, x

    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        print bernoulli_lookup(i)
        s += num / fact * bernoulli(dp.Decimal(i)) * dp.Decimal('-4')**(i / 2) * \
                (1 - dp.Decimal('4')**(i / 2))

    dp.getcontext().prec -= 2
    
    return +s

###
###    miscellaneous mathematical functions.
###

def factorial(n):
    if (n < 1):
        return 1

    else:
        return n * factorial(n - 1)

def combination(m, k):
    if (k <= m):
        return factorial(m)/(factorial(k) * factorial(m - k))

    else:
        return 0

def bernoulli(m):
    if (m == 0):
        return 1

    else:
        t = 0

        for k in range(0, m):
            t += combination(m, k) * bernoulli(k) / (m - k + 1)

        return 1 - t

def bernoulli_lookup(m):

    bernoulli_numerator = [1, 1, 1, 0, -1, 0, 1, 0, -1, 0, 5, 0, -691, 0, 7, 0, -3617, 0, 
        43867, 0, -174611, 0, 854513, 0, -236364091, 0, 8553103, 0, -23749461029, 0, 
        8615841276005, 0, -7709321041217, 0, 2577687858367, 0, -26315271553053477373, 0,
        2929993913841559, 0, -261082718496449122051]

    bernoulli_denominator = [1, 2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 
        798, 1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6, 1, 1919190, 1, 6, 
        1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66, 1, 1590, 1, 798, 1, 870, 1, 354, 
        1, 56786730, 1]


    num = dp.Decimal(bernoulli_numerator[m]) / dp.Decimal(bernoulli_denominator[m])
    return num
