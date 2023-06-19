import sympy as smp

# Define x as a real variable for the function
x = smp.symbols('x', real = True)

# Create the function x^2
function = x ** 2

# Send to integrate the function in terms of x
result = smp.integrate(function, x)

# Pretty print to see it more clearly
smp.pprint(result)

print(smp.latex(result))

