import inspect

# Save a reference to the original built-in print function
original_print = print

def my_print(*args, **kwargs):
    # Get the caller's frame object
    frame = inspect.currentframe().f_back

    # Get the line number of the current line in the caller's frame object
    line_number = frame.f_lineno

    # Call the original print function to output the arguments
    original_print(f"Line {line_number}:", *args, **kwargs)

# Redefine the built-in print function as my_print
print = my_print