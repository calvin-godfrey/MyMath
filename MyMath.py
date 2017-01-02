from __future__ import division #Allows things to be stored as ints but still able to do math
                                #I should probably be using Python 3 for this anyway but too late now.
import re #Used to figure out what class fits a function

def Function(func): #This is sort of a pseudo class, in that it returns a new object
    """
    It takes in a function in basically any form and initializes it as
    a more specialized class. The input is a string in standard mathematical
    form, like it was to be evaluated using normal Python math.
    """
    terms = re.compile("-?(?:(?:(?:\\d+/\\d+|\\d+)\\*)?x(?:\\*\\*-?(?:\\d+/\\d+|\\d+))?|\\d+)")#n*x**m
    terms = terms.findall(func)
    if len(terms) > 1:
        return Polynomial(terms)
    if len(terms) == 1:
        return Term(terms[0])

def simplify(num, denom):
    curr_max = 0
    for i in range(2, int(min(abs(num), abs(denom)))):
        if num%i==0 and denom%i==0:
            curr_max = i
    try:
        num /= curr_max
        denom /= curr_max
    except: #curr_max is still 0
        pass
    return [int(num), int(denom)]

class Term(object):
    """This is the base class for functions of the form n*x**m, can be initialized directly
        or through the Function function above
    It can be initialized in two ways:
    1. A string of the form n*x**m, where n and m are numbers and x can be any letter. (But really only x)
        Note that this is how it would be typed if it were to be evaluated in python.
    2. Keyward arguements, using coefficient, variable, and exponent. Defaults are
        0, x, 0, respectively.
    To make a term a number, with no exponent, either omment the exponent keyward or
        make the exponent 0 with input as a String.
    """
    def __init__(self, *args, **kwargs):
        self.coe_numerator = 0 #These are the default values
        self.coe_denominator = 1
        self.variable = 'x'
        self.exp_numerator = 0
        self.exp_denominator = 1
        if len(args)==0:
            if len(kwargs)==0:
                pass
            else:
                for name, value in kwargs.items():
                    if name == 'coefficient':
                        if isinstance(value, int):
                            self.coe_numerator = value
                            continue
                        value = str(value)
                        if "/" in value:
                            self.coe_numerator, self.coe_denominator = [int(float(i)) for i in value.split("/")]
                        else:
                            self.coe_numerator = int(value)
                    elif name == "variable":
                        self.variable = value
                    elif name == "exponent":
                        if isinstance(value, int):
                            self.exp_numerator = value
                            continue
                        value = str(value)
                        if "/" in value:
                            self.exp_numerator, self.exp_denominator = [int(float(i)) for i in value.split("/")]
                        else:
                            self.exp_numerator = int(value)
        elif len(args) > 1:
            raise TypeError("Term object must be constructed with exactly one string value")
        else:
            if isinstance(args[0], str):
                try:
                    temp_base, exponent = args[0].split('**')
                    if "/" in exponent:
                        self.exp_numerator, self.exp_denominator = [int(i) for i in exponent.split('/')]
                    else:
                        self.exp_numerator = int(exponent)
                        self.exp_denominator = 1
                    self.variable = temp_base[-1]
                    coefficient = temp_base[:-2]
                    if "/" in coefficient:
                        self.coe_numerator, self.coe_denominator = [int(i) for i in coefficient.split("/")]
                    else:
                        try:
                            self.coe_numerator = int(coefficient)
                        except: #It's just x
                            self.coe_numerator, self.coe_denominator = 1, 1
                except ValueError: #No exp (**)
                    value = args[0]
                    if 'x' in value:
                        self.exp_numerator, self.exp_denominator = 1, 1
                        try:
                            coefficient = value.split('x')[0]
                            if coefficient[-1]=='*':
                                coefficient = coefficient[:-1]
                        except IndexError: #There is nothing before "x"
                            self.coe_numerator, self.coe_denominator = 1, 1
                            return None
                    else:
                        self.exp_numerator, self.exp_denominator = 0, 1
                        coefficient = value
                    if "/" in coefficient:
                        self.coe_numerator, self.coe_denominator = [int(i) for i in coefficient.split('/')]
                    else:
                        try:
                            self.coe_numerator = int(coefficient)
                        except ValueError:
                            self.coe_numerator = 0

            else:
                raise TypeError("Unnamed parameter to Term constructor must be a string")
        self.coe_numerator, self.coe_denominator = simplify(self.coe_numerator, self.coe_denominator)
        self.exp_numerator, self.exp_denominator = simplify(self.exp_numerator, self.exp_denominator)
        self.check_sign()

    def __repr__(self):
        """Returns a way to display the term that is basically identical to the string input to __init__"""
        if self.coe_numerator == 0:
            return "0"
        if self.exp_numerator == 0:
            if self.coe_denominator == 1:
                return str(self.coe_numerator)
            return "{}/{}".format(self.coe_numerator, self.coe_denominator)
        elif self.exp_denominator == 1:
            if self.coe_denominator == 1:
                if self.coe_numerator == 1:
                    if self.exp_numerator == 1:
                        return self.variable
                    return "{}**{}".format(self.variable, self.exp_numerator)
                if self.exp_numerator == 1:
                    return "{}*{}".format(self.coe_numerator, self.variable)
                return "{}*{}**{}".format(self.coe_numerator, self.variable, self.exp_numerator)
            if self.coe_denominator == -1:
                return "-{}*{}**{}".format(self.coe_numerator, self.variable, self.exp_numerator)
            return "{}/{}*{}**{}".format(self.coe_numerator, self.coe_denominator, self.variable, self.exp_numerator)
        else:
            if self.coe_numerator == 1:
                return "{}**{}/{}".format(self.variable, self.exp_numerator, self.exp_denominator)
            if self.coe_numerator == -1:
                return "-{}**{}/{}".format(self.variable, self.exp_numerator, self.exp_denominator)
            return "{}/{}*{}**{}/{}".format(self.coe_numerator, self.coe_denominator, self.variable, self.exp_numerator, self.exp_denominator)

    def __add__(self, other):
        """Adds another term to the current one, returning a new one. Only works if the exponents are the same"""
        if not(isinstance(other, Term)):
            raise TypeError("Only a term can be added to another term")
        if self.exp_numerator==other.exp_numerator and self.exp_denominator == other.exp_denominator:
            new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
            new_coe_numerator = new_coe_numerator*other.coe_denominator + other.coe_numerator*new_coe_denominator #Don't have to check type because it has for both of these already
            new_coe_denominator *= other.coe_denominator
            new_coe_numerator, new_coe_denominator = simplify(new_coe_numerator, new_coe_denominator)
            if new_coe_numerator == 0:
                new_exp_numerator = 0
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        else:
            raise ValueError("Exponents of terms must be the same to add. Create a polynomial instead.")

    def __sub__(self, other):
        """Subtract two terms and return a new one"""
        if not(isinstance(other, Term)):
            raise TypeError("Only a term can be subtracted from another term")
        if self.exponent==other.exponent:
            new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
            new_coe_numerator = new_coe_numerator*other.coe_denominator - other.coe_numerator*new_coe_denominator
            new_coe_denominator *= other.coe_denominator
            new_coe_numerator, new_coe_denominator = simplify(new_coe_numerator, new_coe_denominator)
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        else:
            raise ValueError("Exponents of terms must be the same to subtract. Create a polynomial instead.")

    def __mul__(self, other):
        """Return a new term, this times another or a number"""
        new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
        if isinstance(other, Term):
            new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
            new_coe_numerator *= other.coe_numerator
            new_coe_denominator *= other.coe_denominator
            new_coe_numerator, new_coe_denominator = simplify(new_coe_numerator, new_coe_denominator)
            new_exp_numerator = new_exp_numerator*other.exp_denominator + other.exp_numerator*new_exp_denominator
            new_exp_denominator *= other.exp_denominator
            new_exp_numerator, new_exp_denominator = simplify(new_exp_numerator, new_exp_denominator)
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        elif isinstance(other, (int, float)):
            new_coe_numerator *= other
            new_coe_numerator, new_coe_denominator = simplify(new_coe_numerator, new_coe_denominator)
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        else:
            raise TypeError("Terms can only be multiplied by other terms or numbers")

    def __pow__(self, other):
        if isinstance(other, int):
            new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
            new_coe_numerator = new_coe_numerator**other
            new_coe_denominator = new_coe_denominator **other
            new_exp_numerator = new_exp_numerator*other
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))

        raise NotImplementedError #Not yet

    def __truediv__(self, other):
        """Return a new term, this term divided by another or a number"""
        self.check_sign()
        new_exp_numerator, new_exp_denominator, new_coe_numerator, new_coe_denominator = self.exp_numerator, self.exp_denominator, self.coe_numerator, self.coe_denominator
        if isinstance(other, Term):
            other.check_sign()
            if self.exp_numerator==other.exp_numerator and self.exp_denominator==other.exp_denominator and self.coe_denominator == other.coe_denominator and self.coe_numerator == other.coe_numerator:
                return Term("1")
            new_coe_numerator *= other.coe_denominator
            new_coe_denominator *= other.coe_numerator
            new_coe_numerator, new_coe_denominator = simplify(new_coe_numerator, new_coe_denominator)
            new_exp_numerator = new_exp_numerator*other.exp_denominator - other.exp_numerator*new_exp_denominator
            new_exp_denominator *= other.exp_denominator
            new_exp_numerator, new_exp_denominator = simplify(new_exp_numerator, new_exp_denominator)
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        if isinstance(other, (int, float)):
            new_coe_denominator *= other
            return Term(coefficient="{}/{}".format(new_coe_numerator, new_coe_denominator), exponent="{}/{}".format(new_exp_numerator, new_exp_denominator))
        else:
            raise TypeError("Terms can only be divided by other terms or numbers")

    def evaluate(self, num):
        """Return the value of this term at a given value of its variable"""
        return self.coe_numerator/self.coe_denominator * num**(self.exp_numerator/self.exp_denominator)

    def derivative(self):
        """Return a new term that represents the slope of this term for all values of the variable"""
        if self.exp_numerator == 0:
            return Term(coefficient=0)
        temp =  Term(coefficient="{}/{}".format(self.coe_numerator*self.exp_numerator, self.coe_denominator*self.exp_denominator), \
                    variable=self.variable, exponent="{}/{}".format(self.exp_numerator-self.exp_denominator, self.exp_denominator))
        temp.coe_numerator, temp.coe_denominator = simplify(temp.coe_numerator, temp.coe_denominator)
        temp.exp_numerator, temp.exp_denominator = simplify(temp.exp_numerator, temp.exp_denominator)
        return temp

    def indefinite_integral(self):
        """Returns a new term that represents the indefinite integral of this term"""
        if (self.exp_denominator == -1 and self.exp_numerator == 1) or (self.exp_denominator == 1 and self.exp_numerator == -1):
            raise NotImplementedError("Natural log does not exist yet")
        temp = Term(coefficient="{}/{}".format(self.coe_numerator*self.exp_denominator, self.coe_denominator*(self.exp_numerator+self.exp_denominator)),\
                    variable=self.variable, exponent="{}/{}".format(self.exp_numerator+self.exp_denominator, self.exp_denominator))
        temp.coe_numerator, temp.coe_denominator = simplify(temp.coe_numerator, temp.coe_denominator)
        temp.exp_numerator, temp.exp_denominator = simplify(temp.exp_numerator, temp.exp_denominator)
        return temp

    def definite_integral(self, lower, upper):
        """Returns the area under the curve of this term from the lower bound to upper bounds"""
        indefinite = self.indefinite_integral()
        return indefinite.evaluate(upper)-indefinite.evaluate(lower)

    def check_sign(self):
        if self.coe_denominator < 0:
            self.coe_denominator *= -1
            self.coe_numerator *= -1
        if self.exp_denominator < 0:
            self.exp_denominator *= -1
            self.exp_numerator *= -1

class Polynomial(object):
    """This is the class for a series of terms added together"""
    def __init__(self, terms):
        self.terms = []
        if all([isinstance(term, Term) for term in terms]):
            self.terms = terms
            return None
        for term in terms:
            temp_term = Term(term)
            for self_term in self.terms:
                if self_term.exp_numerator == temp_term.exp_numerator and self_term.exp_denominator == temp_term.exp_denominator:
                    self_term += temp_term
                    break
            else: #Term with that exponent is not already in the list
                self.terms.append(temp_term)
        for term in self.terms:
            term.check_sign()
        self.terms = sorted(self.terms, key=lambda i:i.exp_numerator/i.exp_denominator)[::-1]

    def __repr__(self):
        self.combine()
        return ''.join([repr(term) if repr(term)[0]=="-" or index==0 else "+"+repr(term) for index,term in enumerate(self.terms)])

    def __add__(self, other):
        new_arr = self.terms
        if isinstance(other, (Term, int)):
            if isinstance(other, int):
                other = Term(coefficient=other)
            for term in new_arr:
                if term.exp_numerator/term.exp_denominator == other.exp_numerator/other.exp_denominator:
                    term += other
                    break
            else: #The polynomial doesn't contain a term with a matching power
                new_arr.append(other)
        elif isinstance(other, Polynomial):
            for other_term in other.terms:
                for this_term in new_arr: #TODO: make this more efficient somehow
                    if other_term.exp_numerator/other_term.exp_denominator == this_term.exp_numerator/this_term.exp_denominator:
                        this_term += other_term
                        break
                new_arr.append(other_term)
        else:
            raise TypeError("Only a polynomial, term, or int can be added to a polynomial")
        for term in new_arr:
            term.check_sign()
        return Polynomial(new_arr)

    def __sub__(self, other):
        neg_other_arr = []
        curr_terms = Polynomial(self.terms)
        for term in other.terms:
            term *= -1
            neg_other_arr.append(term)
        neg_other = Polynomial([repr(i) for i in neg_other_arr])
        curr_terms = curr_terms+neg_other
        curr_terms.combine()
        return Polynomial(curr_terms.terms)

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            new_arr = []
            for other_term in other.terms:
                for this_term in self.terms:
                    new_arr.append(other_term*this_term)
            return Polynomial(new_arr)
        if isinstance(other, (Term, int)):
            new_arr = self.terms
            if isinstance(other, int):
                other = Term(coefficient=other)
            for term in new_arr:
                term *= other
            return Polynomial(new_arr)
        else:
            raise TypeError("Polynomial can only be multiplied by int, term, or polynomial")

    def __truediv__(self, other): #THIS IS JUST COMPLETELY WRONG
        if isinstance(other, Polynomial): #Assuming there is no remainder in division
            copy_self = Polynomial(self.terms)
            ans = []
            while len(copy_self.terms)>1:
                div = copy_self.terms[0]/other.terms[0]
                mult_other = Polynomial([term*div for term in other.terms])
                copy_self -= mult_other
                copy_self.combine()
                copy_self.arrange()
                ans.append(div)
            ans = Polynomial(ans)
            ans.arrange()
            ans.combine()
            ans_str = ''.join([repr(i) if repr(i)[0]=='-' or index==0 else "+"+repr(i) for index, i in enumerate(ans.terms)])
            factor = re.compile("-?(?:(?:(?:\\d+/\\d+|\\d+)\\*)?x(?:\\*\\*-?(?:\\d+/\\d+|\\d+))?|\\d+)")
            search = factor.findall(ans_str)
            return Polynomial(ans.terms)

        elif isinstance(other, (Term, int)): #Assuming there is no remainder in division
            new_arr = list(self.terms)
            if isinstance(other, int):
                other = Term(coefficient="{}/1".format(other))
            for term in new_arr:
                term /= other
            return Polynomial(new_arr)
        else:
            raise TypeError("Polynomial can only be divided by int, term, or polynomial")

    def __pow__(self, other):
        if isinstance(other, int):
            if other > 0:
                temp = Polynomial(self.terms)
                ans = Polynomial(self.terms)
                for i in range(1, other):
                    ans = ans*temp
                return ans
        raise NotImplementedError

    def evaluate(self, num):
        return sum([term.evaluate(num) for term in self.terms])

    def derivative(self):
        return Polynomial([term.derivative() for term in self.terms])

    def indefinite_integral(self):
        return Polynomial([term.indefinite_integral() for term in self.terms])

    def definite_integral(self, lower, upper):
        return sum([term.definite_integral(lower, upper) for term in self.terms])

    def combine(self):
        for i in xrange(len(self.terms)):
            for j in xrange(i+1, len(self.terms)):
                terms1 = self.terms[i]
                try:
                    terms2 = self.terms[j]
                except IndexError:
                    self.terms = filter(lambda x:repr(x)!='0', self.terms)
                    return self
                try:
                    terms1 = terms1+terms2
                    self.terms[i] = terms1
                    self.terms.pop(j)
                    j -= 1
                except ValueError:
                    continue
        self.terms = filter(lambda x:repr(x)!='0', self.terms)
        return self

    def factor(self):
        min_exp_numerator = 1e10
        min_exp_denominator = 1
        new_arr = list(self.terms)
        final = []
        ans = ""
        for term in self.terms:
            if term.exp_numerator/term.exp_denominator < min_exp_numerator/min_exp_denominator:
                min_exp_numerator = term.exp_numerator
                min_exp_denominator = term.exp_denominator
        if min_exp_numerator < 1e10:
            alone = Term("x**{}/{}".format(min_exp_numerator, min_exp_denominator))
            for term in new_arr:
                term = term/Term("x**{}/{}".format(min_exp_numerator, min_exp_denominator))
                final.append(term)
            if repr(alone)=='1':
                alone = ""
        else:
            alone = ""
        if alone=="":
            pass
        else:
            ans += repr(alone)
        copy_self = Polynomial(final)
        for i in range(-100, 101): #Reasonable range for factors
            try:
                if copy_self.evaluate(i) == 0:
                    if i<0:
                        temp = Function("x+"+str(abs(i)))
                    else:
                        temp = Function("x-"+str(i))
                    copy_self = copy_self/temp
                    if len(ans):
                        ans += "*({})".format(repr(temp))
                    else:
                        ans += "({})".format(repr(temp))
                    if repr(copy_self)=='1' or len(copy_self.terms) < 1:
                        break
            except AttributeError:
                break
        return ans

    def arrange(self):
        self.terms = sorted(self.terms, key=lambda i:i.exp_numerator/i.exp_denominator)[::-1]
        return self

    def check_sign(self):
        self.terms = [term.check_sign() for term in self.terms]
        return self

    def __eq__(self, other):
        self.check_sign()
        self.arrange()
        other.check_sign()
        other.arrange()
        return self.exp_numerator==other.exp_numerator and self.exp_denominator==other.exp_denominator and self.coe_numerator==other.exp_numerator and self.coe_denominator==other.coe_denominator


class Rational(object):
    """Class for functions that have x's in numerator and denominator"""
    def __init__(self, numerator, denominator):
        if isinstance(numerator, (Polynomial, Term)):
            self.numerator = numerator
        else:
            self.numerator = Function(numerator)
        if isinstance(denominator, (Polynomial, Term)):
            self.denominator = denominator
        else:
            self.denominator = Function(denominator)

b = Function("x**2+9*x-10")
e = Function("x**3+11*x**2+8*x-20")
d = Function("x**2+2*x")
c = Function("x+2")
print e.factor()
print b.factor()
print d.factor()
print d.factor()
print d/c
