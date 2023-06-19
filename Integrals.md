# **Integral examples**

Examples of how to write functions and integrate them in `integrals.py` using `sympy`

To obtain the `Latex` for the function or result you can execute

```
print(smp.latex(function))
```

Where function can be passed as the result of the integration or the function itself.

This can help view the result better if the `pprint()` function is too complex

### **First example**

Calculate the integral of $x^2$

```
# Define x as a real variable for the function
x = smp.symbols('x', real = True)

# Create the function x^2
function = x ** 2

# Send to integrate the function in terms of x
result = smp.integrate(function, x)
```
Result: $\frac{x^{3}}{3}$
