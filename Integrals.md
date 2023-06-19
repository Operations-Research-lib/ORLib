# **Integral examples**

Examples of how to write functions and integrate them in `integrals.py` using `sympy`

To obtain the `Latex` for the function or result you can execute

```
print(smp.latex(function))
```

Where function can be passed as the result of the integration or the function itself.

This can help view the result better if the `pprint()` function is too complex

### **Basic sympy notation**

1. $x^n$ = `x ** n`

3. $\log_{y}(x)$ = `smp.log(x, y)`

4. $e^x$ = `smp.exp(x)`

2. $\sqrt{x}$ = `smp.sqrt(x)`

5. $\frac{a}{b}$ = `smp.Rational(a,b)`

6. $\cos(x)$ = `smp.cos(x)`

7. $\sin(x)$ = `smp.sin(x)`

8. $\tan(x)$ = `smp.tan(x)`


### **Notes:**

+ $x^n$ does not require `smp` beforehand as it is base python notation

+ smp.Rational(a,b) is used for the number to be seen as a fraction instead of a floating point number

+ The rest of the trigonometric functions can be obtained by using: 

```
smp.asin(x) 
smp.acos(x) 
smp.atan(x)
smp.acot(x)
smp.sinh(x)
smp.cosh(x)
smp.tanh(x)
smp.coth(x)
smp.asinh(x)
smp.acosh(x)
smp.atanh(x)
smp.acoth(x)
```

## **Examples**

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

### **Second example**

Calculate the integral of $\frac{x}{2}$

```
# Define x as a real variable for the function
x = smp.symbols('x', real = True)

# Create the function
function = x / 2

# Send to integrate the function in terms of x
result = smp.integrate(function, x)
```

Result: $\frac{x^{2}}{4}$

### **Third example**

Calculate the integral of $\frac{1}{\sqrt{x} \left(x + 1\right)}$

```
# Define x as a real variable for the function
x = smp.symbols('x', real = True)

# Create the function
function = 1 / (smp.sqrt(x) * (x + 1))

# Send to integrate the function in terms of x
result = smp.integrate(function, x)
```

Result: $2 \operatorname{atan}{\left(\sqrt{x} \right)}$

