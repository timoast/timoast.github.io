---
title: Using python decorators
author: Tim Stuart
date: '2017-04-09'
comments: true
---

Yesterday I wrote my first python decorator. Decorators have always seemed a bit mysterious to me, but having finally written one I can see a bit better how they work. This is the decorator I wrote:

<!--break-->

```python
def log_info(func):
    def wrapper(args):
        print("Function {} called with the following arguments:\n".format(func.__name__))
        for arg in vars(args):
            print(str(arg) + '\t' + str(getattr(args, arg)))
        t1 = time.time()
        func(args)
        t2 = time.time()
        elapsed = [round(x, 2) for x in divmod(t2-t1, 60)]
        print("\nFunction completed in  {} m {} s\n".format(elapsed[0], elapsed[1]))
    return wrapper
```

It just takes a function and returns a new function that is a 'wrapped' version of the original function. This wrapped version and prints the name of the function, the list of arguments the function received, runs the original function, then prints how long it took to complete.

A decorator is a function that takes a function as its argument and returns a function. A decorator can be used to alter what a function does without having to change the code of the function itself. We can see that in my example above, the function `log_info` takes a function and defines a new function within its scope, and then returns this new function.

## How decorators work

This is a lot to take in, so let's start at the beginning.

### What is a function

A functions is something generates a value based on given arguments. A function is an object. Arguments to a function are also objects. A function can accept a function as its argument.

### Nested functions

As seen in our above example, functions can be defined inside another function. This is known as a nested function.


```python
def outer():
    def inner():
        return 1
    return inner

nested_example = outer()
print(nested_example())
```

```
## 1
```

### What is a closure

A closure can be produced by a nested function. Closures remember their enclosing scope from when they were defined. In this way, closures can have values hard-coded determined by their enclosing scope at the time they were defined. This can avoid the use of global variables. For example:


```python
def outer(x):
    def inner():
        return x
    return inner

closure1 = outer(1)
closure2 = outer(2)
print(closure1())
print(closure2())
```

```
## 1
## 2
```

`closure1` returns 1, even though no arguments are given. This is because the value as set when the function was assigned using a closure.

### What is a decorator

Now, back to our original question. A decorator is a function that takes and returns another function. They use nested functions and closures. Take this as a minimal example:


```python
def my_decorator(func):
    def inner():
        return func()
    return inner
```

The decorator make use of the fact that the inner function has access to objects in the enclosing scope. That's why we don't need to pass `inner()` the argument `func`. Now, if we were to use this function as a decorator, we would do this:


```python
my_function = my_decorator(my_function)
print(my_function())
```

```
## 1
```

`my_function` has now been 'decorated' with `my_decorator`. That is, it's been re-assigned as the output of `my_decorator`. In our example, the decorator doensn't actually do anything. Here's a more useful example:


```python
def multiply(func):
    def inner():
        return func() * 2
    return inner

my_function = multiply(my_function)
print(my_function())
```

```
## 2
```

Now `my_function` returns 2 not 1. If we wanted to be able to specify variables in the decorator function, we do so using another enclosing function (basically, we make a decorator that returns a decorator):


```python
def multiply(multiplier):
    def outer(func):
        def inner():
            return func() * multiplier
        return inner
    return outer

my_function = multiply(3)(my_function)
print(my_function())
```

```
## 3
```

This syntax (`multiply(3)(my_function)`) does look a bit strange. Python has another way to decorate functions, using the `@decorator` syntax before the function definition. If we take our first example, with no parameters, would could have written it like this: 


```python
@my_decorator
def my_function():
    return 1

print(my_function())
```

```
## 1
```

And if we look at our second example, with parameters:


```python
@multiply(3)
def my_function():
    return 1

print(my_function())
```

```
## 3
```

Looks much better!

Let's go back to my original example, the decorator I wrote yesterday. Hopefully now it makes a bit more sense. It takes a function, then returns a function that runs the input function plus a bit more stuff. We can use it to decorate whatever functions we want. Below is an example. I've just changed how it prints out the function arguments, as the way I wrote it originally deals with command-line options:


```python
import time

def log_info(func):
    def wrapper(args):
        print("Function {} called with the following arguments:\n{}\n".format(func.__name__, args))
        t1 = time.time()
        output = func(args)
        t2 = time.time()
        elapsed = [round(x, 2) for x in divmod(t2-t1, 60)]
        print("Function completed in  {} m {} s\n".format(elapsed[0], elapsed[1]))
        return output
    return wrapper

@log_info
def my_function(number):
    return number * 3

x = my_function(2)
print(x)
```

```
## Function my_function called with the following arguments:
## 2
## 
## Function completed in  0.0 m 0.0 s
## 
## 6
```
