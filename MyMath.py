from __future__ import division #Allows things to be stored as ints but still able to do math
                                #I should probably be using Python 3 for this anyway but too late now.
import re #Used to figure out what class fits a function

def Function(func): #This is sort of a pseudo class, in that it returns a new object
    """
    It takes in a function in basically any form and initializes it as
    a more specialized class. The input is a string in standard mathematical
    form, like it was to be evaluated using normal Python math.
    """
    terms = re.compile("-?(?:\\d+/\\d+|\\d+)\\*x\\*\\*-?(?:\\d+/\\d+|\\d+)")#n*x**m
    terms =  terms.findall(func)
    if len(terms) > 1:
        return Polynomial(terms)
    if len(terms) == 1:
        return Term(func)

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
                            self.coe_numerator, self.coe_denominator = [int(i) for i in value.split("/")]
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
                            self.exp_numerator, self.exp_denominator = [int(i) for i in value.split("/")]
                        else:
                            self.exp_numerator = int(value)
        elif len(args) > 1:
            raise TypeError("Term object must be constructed with exactly one string value")
        else:
            if isinstance(args[0], str):
                temp_base, self.exponent = args[0].split('**')
                if "/" in self.exponent:
                    self.exp_numerator, self.exp_denominator = [int(i) for i in self.exponent.split('/')]
                else:
                    self.exp_numerator = int(self.exponent)
                    self.exp_denominator = 1
                self.variable = temp_base[-1]
                self.coefficient = temp_base[:-2]
                if "/" in self.coefficient:
                    self.coe_numerator, self.coe_denominator = [int(i) for i in self.coefficient.split("/")]
                else:
                    self.coe_numerator = int(self.coefficient)
            else:
                raise TypeError("Unnamed parameter to Term constructor must be a string")
        self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
        self.exp_numerator, self.exp_denominator = self.reduce(self.exp_numerator, self.exp_denominator)

    def __repr__(self):
        """Returns a way to display the term that is basically identical to the string input to __init__"""
        if self.coe_numerator==0:
            return Term(coefficient="0")
        if self.exp_numerator == 0:
            if self.coe_numerator == 1:
                return str(self.coe_numerator)
            return "{}/{}".format(self.coe_numerator, self.coe_denominator)
        elif self.exp_denominator == 1:
            if self.coe_denominator == 1:
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
        """Adds another term to the current one, only works if the exponents are the same"""
        if not(isinstance(other, Term)):
            raise TypeError("Only a term can be added to another term")
        if self.exponent==other.exponent:
            self.coe_numerator = self.coe_numerator*other.coe_denominator + other.coe_numerator*self.coe_denominator #Don't have to check type because it has for both of these already
            self.coe_denominator *= other.coe_denominator
            self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            return self
        else:
            raise ValueError("Exponents of terms must be the same to add. Create a polynomial instead.")

    def __sub__(self, other):
        """Subtract two terms and change this one"""
        if not(isinstance(other, Term)):
            raise TypeError("Only a term can be subtracted from another term")
        if self.exponent==other.exponent:
            self.coe_numerator = self.coe_numerator*other.coe_denominator - other.coe_numerator*self.coe_denominator
            self.coe_denominator *= other.coe_denominator
            self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            return self
        else:
            raise ValueError("Exponents of terms must be the same to subtract. Create a polynomial instead.")

    def __mul__(self, other):
        """Multiply a Term by different things, either term or number"""
        if isinstance(other, Term):
            self.coe_numerator *= other.coe_numerator
            self.coe_denominator *= other.coe_denominator
            self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            self.exp_numerator = self.exp_numerator*other.exp_denominator + other.exp_numerator*self.exp_denominator
            self.exp_denominator *= other.exp_denominator
            self.exp_numerator, self.exp_denominator = self.reduce(self.exp_numerator, self.exp_denominator)
            return self
        elif isinstance(other, (int, float)):
            self.coe_numerator *= other
            self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            return self
        else:
            raise TypeError("Terms can only be multiplied by other terms or numbers")

    def __pow__(self, other):
        raise NotImplementedError #Not yet

    def __truediv__(self, other):
        """Divide a term by other things, either term or number"""
        if isinstance(other, Term):
            self.coe_numerator *= other.coe_denominator
            self.coe_denominator *= other.coe_numerator
            self.coe_numerator, self.coe_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            self.exp_numerator = self.exp_numerator*other.exp_denominator - other.exp_numerator*other.exp_numerator
            self.exp_denominator *= other.exp_denominator
            self.exp_numerator, self.exp_denominator = self.reduce(self.coe_numerator, self.coe_denominator)
            return self
        if isinstance(other, (int, float)):
            self.coe_denominator *= other
            return self
        else:
            raise TypeError("Terms can only be divided by other terms or numbers")

    def reduce(self, num, denom):
        curr_max = 0
        for i in range(2, int(min(abs(num), abs(denom)))):
            if num%i==0 and denom%i==0:
                curr_max = i
        try:
            num /= curr_max
            denom /= curr_max
        except:
            pass
        return [int(num), int(denom)]

    def evaluate(self, num):
        """Return the value of this term at a given value of its variable"""
        return self.coe_numerator/self.coe_denominator * num**(self.exp_numerator/self.exp_denominator)

    def derivative(self):
        """Return a new term that represents the slope of this term for all values of the variable"""
        if self.exp_numerator == 0:
            return Term(coefficient=0)
        temp =  Term(coefficient="{}/{}".format(self.coe_numerator*self.exp_numerator, self.coe_denominator*self.exp_denominator), \
                    variable=self.variable, exponent="{}/{}".format(self.exp_numerator-self.exp_denominator, self.exp_denominator))
        temp.coe_numerator, temp.coe_denominator = temp.reduce(temp.coe_numerator, temp.coe_denominator)
        temp.exp_numerator, temp.exp_denominator = temp.reduce(temp.exp_numerator, temp.exp_denominator)
        return temp

    def indefinite_integral(self):
        """Returns a new term that represents the indefinite integral of this term"""
        if self.exponent == -1:
            raise NotImplementedError("Natural log does not exist yet")
        temp = Term(coefficient="{}/{}".format(self.coe_numerator*self.exp_denominator, self.coe_denominator*(self.exp_numerator+self.exp_denominator)),\
                    variable=self.variable, exponent="{}/{}".format(self.exp_numerator+self.exp_denominator, self.exp_denominator))
        temp.coe_numerator, temp.coe_denominator = temp.reduce(temp.coe_numerator, temp.coe_denominator)
        temp.exp_numerator, temp.exp_denominator = temp.reduce(temp.exp_numerator, temp.exp_denominator)
        return temp

    def definite_integral(self, lower, upper):
        """Returns the area under the curve of this term from the lower bound to upper bounds"""
        indefinite = self.indefinite_integral()
        return indefinite.evaluate(upper)-indefinite.evaluate(lower)

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
                if self_term.exponent == temp_term.exponent:
                    self_term += temp_term
                    break
            else:
                self.terms.append(temp_term)
        self.terms = sorted(self.terms, key=lambda i:i.exponent)

    def __repr__(self):
        return ''.join([repr(term) if repr(term)[0]=="-" or index==0 else "+"+repr(term) for index,term in enumerate(self.terms)])

    def __add__(self, other):
        if isinstance(other, (Term, int)):
            if isinstance(other, int):
                other = Term(coefficient=other)
            for term in self.terms:
                if term.exponent == other.exponent:
                    term += other
                    break
            else: #The polynomial doesn't contain a term with a matching power
                self.terms.append(other)
        elif isinstance(other, Polynomial):
            for other_term in other.terms:
                for this_term in self.terms: #TODO: make this more efficient somehow
                    if other_term.exponent == this_term.exponent:
                        this_term += other_term
                        break
                self.terms.append(other_term)
        else:
            raise TypeError("Only a polynomial, term, or int can be added to a polynomial")
        return self

    def __sub__(self, other):
        neg_other_arr = []
        for term in other.terms:
            term *= -1
            neg_other_arr.append(term)
        neg_other = Polynomial(''.join(neg_other_arr))
        self += neg_other
        return self

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            new_arr = []
            for other_term in other.terms:
                for this_term in self.terms:
                    new_arr.append(other_term*this_term)
            self = Polynomial(''.join(new_arr))
            return self
        if isinstance(other, (Term, int)):
            if isinstance(other, int):
                other = Term(coefficient=other)
            for term in self.terms:
                term *= other
            return self
        else:
            raise TypeError("Polynomial can only be multiplied by int, term, or polynomial")

    def __truediv__(self, other):
        if isinstance(other, Polynomial):
            new_arr = []
            for other_term in other.terms:
                for this_term in self.terms:
                    new_arr.append(this_term/other_term)
            self = Polynomial(''.join(new_arr))
            return self
        elif isinstance(other, (Term, int)):
            if isinstance(other, int):
                other = Term(coefficient=other)
            for term in self.terms:
                term /= other
            return self
        else:
            raise TypeError("Polynomial can only be divided by int, term, or polynomial")

    def evaluate(self, num):
        return sum([term.evaluate(num) for term in self.terms])

    def derivative(self):
        return Polynomial([term.derivative() for term in self.terms])

    def indefinite_integral(self):
        return Polynomial([term.indefinite_integral() for term in self.terms])

    def definite_integral(self, lower, upper):
        return sum([term.definite_integral(lower, upper) for term in self.terms])

'''x = Function("4*x**2")
print x.evaluate(6)
print x.derivative()
print x.indefinite_integral()
print x.definite_integral(0, 10)
y = Term(coefficient=-5, exponent="2") #Different way to construct it, not all kwargs required, can be string or int or float
x *= y #Right now it pretends that all theh variables are the same and ignores them in math
print x'''

z = Function("-8/3*x**3/2+5/3*x**7/8")
print z.derivative()
print z.indefinite_integral()
