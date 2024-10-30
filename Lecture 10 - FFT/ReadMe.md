This code can be used to price a barrier option with discrete monitoring under any Lévy, exactly like we previously did with MC. And it's a general method : the only thing we have to change if we want to use another kind of Lévy is the computation of the characteristic function. Nothing else. \\

**Warning :** it's important to understand why we shift : MATLAB wants 0->2N, we want -N->N-1.
