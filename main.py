import sympy as sp

class Numerical:

    def __init__(self):
        # unequal interval
        self.x = [5, 6, 9, 11]
        self.fx = [12, 13, 14, 16]

        # equal interval
        self.x2 = [5, 10, 15, 20]
        self.fx2 = [90, 105, 135, 180]

    def lagranges_interpolation(self):
        # unequal interval
        t = sp.symbols("t")
        y = 0  

        for i in range(len(self.x)):
            L_i = 1  
            for j in range(len(self.x)):
                if i != j:
                    L_i *= (t - self.x[j]) / (self.x[i] - self.x[j])

            y += L_i * self.fx[i]

        return sp.simplify(y)
    
    def newton_divided_difference(self):
        t = sp.symbols("t")
        
        n = len(self.x)
        divided_diff = [[0] * n for _ in range(n)] # for all the dell x list storage
        
        for i in range(n):
            divided_diff[i][0] = self.fx[i]
        
        
        for j in range(1, n):  # j  column index
            for i in range(n - j):  # i row index
                # it works becs it takes j-1 column and we are progressing towards divided diff
                divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (self.x[i + j] - self.x[i])
        # print(divided_diff)
        y = self.fx[0]
        term = 1  # Keeps track of (t - x0)(t - x1)...(t - xi)
        for j in range(1, n):
            term *= (t - self.x[j - 1])
            y += term * divided_diff[0][j]
        
        return sp.simplify(y)
    
    def nfb(self, method="forward"): # newton forward backward helper
        # f(a + h*u) = f(a) + u*del(f(a)) + (u*(u-1)/2)*del2(f(a)) + so on
        u = sp.symbols("u")
        h = self.x2[1] - self.x2[0] # difference
        n = len(self.x2)
        delta = [[0] * n for _ in range(n)] # for all the dell x list storage
        
        for i in range(n):
            delta[i][0] = self.fx2[i]
        
        for j in range(1, n):  # j  column index
            for i in range(n - j):  # i row index
                # it works becs it takes j-1 column and we are progressing towards further delta
                delta[i][j] = (delta[i + 1][j - 1] - delta[i][j - 1])

        if method == "forward":
            polynomial = delta[0][0] 
            u_term = 1  
            for j in range(1, n):  # Add terms one by one
                u_term *= (u - (j - 1)) / j
                polynomial += u_term * delta[0][j]
            return polynomial
        
        # now lets say for the example above if asked the value of polynomial at x =8
        # 8= 5[x0] + (0.6)*5[common difference] . so in the polynomial if u = 0.6 ans can be calculated 
        # how to calculate u  based on numerical input in case of newton forward
        # (input - first_term(x@0))/gap newton forward
        # (input - last_term(x@-1))/gap newton backward
        else:
            polynomial = delta[0][n-1] #start backwards 
            u_term = 1  
            for j in range(1, n):  # Add terms one by one
                u_term *= (u + (j - 1)) / j
                polynomial += u_term * delta[n-j-1][j]
            return polynomial
    
    def newton_forward_interpolation(self, input):
        u = sp.symbols("u")
        eq = self.nfb(method="forward")

        # (input - first_term(x0))/gap newton forward
        value = eq.subs(u, (input-self.x2[0])/(self.x2[1]-self.x2[0])).evalf()

        return value
    
    def newton_backward_interpolation(self, input):
        u = sp.symbols("u")
        eq = self.nfb(method="backward")

        # (input - last_term(x@-1))/gap newton backward
        value = eq.subs(u, (input-self.x2[-1])/(self.x2[1]-self.x2[0])).evalf()

        return value
    
    def newton_forward_diff(self, nth, input):
        #nth is the n th derivative 
        u = sp.symbols("u")
        eq = self.nfb(method="forward")

        derivative = sp.diff(eq, u, nth)

        diff = (self.x2[1]-self.x2[0])

        value = derivative.subs(u, (input-self.x2[0])/diff).evalf()/(diff**nth)

        return value
    
    def newton_backward_diff(self, nth, input):
        #nth is the n th derivative 
        u = sp.symbols("u")
        eq = self.nfb(method="backward")

        derivative = sp.diff(eq, u, nth)

        diff = (self.x2[1]-self.x2[0])

        value = (derivative.subs(u, (input-self.x2[-1])/diff).evalf())/(diff)**nth

        return value


        
print(Numerical().newton_backward_diff(1, 7))
