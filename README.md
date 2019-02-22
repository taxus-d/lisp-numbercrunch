# lisp-numbercrunch

### Contents
currently includes rather naive recursive fft implementaion that
may be useful to learn which directives allow sbcl to heavily optimize code
### Problems
* NO vectorization. Fortran compilers can use it almost out-of-the-box. 
Things are much more complicated with sbcl 
