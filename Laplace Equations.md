*Numerical Methods* 

**D.**   
13 

**D.1.**   
**Parallel** Simulation **of the Time** Independent **Diffusion** Equation 

**Introduction** 

In many cases one is only interested in the final steady state concentration field and not so much in the transient behavior, i.e. the route towards the steady state is not relevant. This may be so because the diffusion is a very fast process in comparison with other processes in the system, and then one may neglect the transient behavior (an example is the diffusion limited aggregation that is presented as a case study at the end of this chapter). 

Setting all time derivatives to zero in the diffusion equation (equation \[4\]) we find the time independent diffusion equation: 

V2c=0.   
\[14\] 

This is the famous *Laplace equation* that is encountered in many places in science and engineering, and is therefore very relevant to study in some more detail. 

Again consider the two-dimensional domain and the discretisation of the previous section. The concentration c no longer depends on the time variable, so we denote the concentration on the grid point here as *cm*. Substituting the finite difference formulation for the spatial derivatives (Equation \[10\]) into Equation \[14\] results in 

*€1,m* \= 1 | C1\+1,m \+*91-1*,m *\+ C1,*m\+1 \+*C1,*m−1\], V(l,m) .   
·   
\[14\] 

Equation \[14\] contains a set of (*L*\+1)x(*L*\+1\) linear equations with (L+1)×(*L*\+1) unknowns. The unknowns are the value of the concentration of all grid points. It is our task to solve this set of equations. These equation can be reformulated as 

**Ax \= b**,   
\[15\] 

\=   
where **A** is a (NXN) matrix and **x** and **b** are *(N*×1) vectors. As an example, let us derive an explicit formulation for **A** in the one-dimensional case. Consider a one\-dimensional domain 0 ≤ x ≤ 1 and a discretisation with *Ex* 1/*L,* such that *c1* \= *c(l\&x*) with 0 ≤1 ≤ *L.* Furthermore, assume that the concentration on the boundaries of the domain are fixed, i.e. *co* \= CO and *c* \= *CL.* The one-dimensional Laplace equation reads a2c*/*dx2 \= 0, and therefore the finite difference equations become   
*Co* 

*C1-1 −2c1 \+C1*\+1 \= 0,   
0≤1≤L. 

In other words 

2 \-1 0   
CO-2c1 \+*c2* \= 0 *C1*\-*2c2*\+*c3*\=0   
\-1 2   
от 

0 \-1 2 \-1 

or 

: 

CL\-3-2CL*\-*2 \+*CL*\-1\=0 

*CL*\-2-2c1-1\+*CL*\=0   
:   
0   
C1   
CO 

C2   
0 

: 

\-1 2 \-1 *CL*\-2   
0 

0   
\-1 2   
CL\-1   
*CL* 

In this example we therefore find that **A** is a tri-diagonal matrix with *N \= L* − 1\. In higher dimensions **A** is still banded but has a somewhat more complicated structure. You may try yourself to find **A** for the two dimensional domain and boundary conditions used in the previous section.   
*Numerical Methods* 

**D.2**.   
14 

Numerically solving sets of linear equations is a huge business and we can only give some important results and discuss some of the implications with respect to parallelism. A more in depth discussion can be found in the appendix, section E, and references therein. 

\=   
One may solve the set of equation **Ax** methods one solves the equations directly within some numerical error. In iterative methods one seeks an approximate solution via some iterative procedure. Table 1 shows the computational complexity for direct and iterative methods, both for banded and dense matrices. For iterative methods it is assumed that the number of iterations needed to find a solution within certain accuracy is much smaller than *N.* Direct methods allow finding the solution within numerical accuracy at the expense of a higher computational complexity. Especially for very large matrix dimensions, as are routinely encountered in today's large scale simulations this large computational complexity becomes prohibitive, even when executed on high-end massively parallel computers. That is the most important reason to consider iterative methods. They are cheaper, typically easier to parallelise, however at the expense only finding approximate solution. Sometimes convergence can be very slow, and iterative methods then require very specialized (preconditioning) techniques to recover convergence. These issues will not be discussed here. We will first give a short survey of some of the direct methods to solve a set of equations, and then turn our attention to iterative methods, applied to the case of the discretised two-dimensional Laplace equation. Much more in-depth discussion can be found in the appendix, section E and reference therein.   
**b** using *direct* or *iterative* methods. In direct 

dense matrix   
banded matrix 

direct method 

iterative method   
O(*N3*) O(N2)   
N 

O(N2) 

O(N) 

*Table 1: Computational complexity for direct and iterative methods*. 

**Direct Methods** 

Gaussian Eliminiation 

The standard way to solve a system of equations is by *Gaussian Elimination.* You may know this method also under the name *matrix sweeping*. The idea is to bring the system of equations in *upper triangular form.* This means that all elements of the matrix below the diagonal are zero. So, the original system of equations 

a1N   
*b1*   
\[h\! 1\\?   
*a NN* 

is transformed to the system   
*N*   
*N* 

a11   
a12   
...   
αν 

0   
a22 

0   
0   
a'N − 1, *N a'NN*   
x1 

X   
N   
*b'*   
*N* 

This upper triangular form is easily solved through *back substitution.* First calculate   
*Numerical Methods* 

XN 

followed by 

*\**N-1   
a'NN 

*BN*\-1-N\-1,N*\*N* 

*"N*\-1,*N*\-1 

etc. In a single expression we find   
15 

*N* 

\=*i*\+1' 

*xi*   
*i* \= N,···,1. 

Our main concern is to transform the original system **Ax** 

form **A'x** matrix   
\=   
**b** into the upper triangular **b'**. This is achieved using Gaussian Eliminiation. Consider the composed   
\= 

a1N   
a11   
*b1* :   
**(A,b)** \= 

a *NN*   
bN   
a N1 

The recipe for Gaussian Elimination is 

a) Determine an element *a,*1 0, proceed with (b). If such element does not exist, **A** is 

singular, stop. 

b) Interchange row *r* and row 1 of **(A, b)**, the result is **(A,** b) . 

c) For *i* \= 2, 3, *N*, subtract multiple 711 \= *ā¡/*ā11 of row 1 from row *i.* This results is 

matrix 

*a* 

(1)   
0   
a (a)TM   
**(A**,**b**) \= 

**A**   
**b** 

0   
*αγι*   
1 bo   
*NN* 

Next, the same procedure is applied to **(A,b**) and repeated, i.e 

**(A,b**) :\= **(A**,b)→**(A®**,b®)→   
(*N*\-1)   
··· → **(A b** (N-1)) \= (U, **c)**.   
, 

The system (**U,c)** has the desired upper triangular form and can be solved easily by back substitution. 

\=   
The element *a,* is called the *pivot* element, and the choice of the pivot is very important. In general it is not a good idea to take e.g. *ar* \= a11. This may lead to loss of accuracy, or even instability of the Gaussian Elimination. This happens is a11 is relatively small compared to the other elements. In that way the multipliers *l11* may become very large and in the subtraction of rows we may end up subtracting large numbers from each other, leading to enormous loss of precision. To avoid such problems we must always apply *partial pivoting*, i.e. choose |*a*,|\= max|a|, i.e. pick the maximum element (in absolute   
*i* 

measure) as the pivot. In general this procedure leads to stable and accurate solutions. For more information on Gaussian Elimination, read section **Error\! Reference source not found**. of the appendix. 

LU Factorization 

Let us take a closer look at the sweeping of the matrix. Note that step (b), the interchanging of rows, can be formulated as   
*Numerical Methods*   
16 

where **P1**   
**(A,b)** \= P1 **(A, b)** 

is a permutation matrix. Consider for example the case where *N* \=4 and we want to exchange row 1 and 2\. In that case we find 

0 1 0 0 

1   
0 0 0   
**P1**   
0 0   
1 0 

0 0   
0 1 

Note that P1P1 \= I, where **I** is the identity matrix. From this we can conclude, and this is   
1 

\-1   
generally true for the permutation matrix, that **P** \= P1, i.e. P1 is its own inverse. Also the actual sweeping in step (c) can be formulated as 

where   
**(A,**b®)\=G1**(A, b**) 

1 0 

\-21 1   
**G1**   
: 

\-INI 0 1 

Note that the inverse of **G1** is immediately found to be 

1   
0   
Ο 

*721*   
1 

INI   
0 1 

1   
Combining this results in (**A**,b) \= **G1P1 (A,**b), which is easily generalized to (A®,b®)\=G ̧P;(A), b)), where P; is a permutation matrix and **G**, is 

1 

1 

: 

0   
\-INJ 

With all these definitions we can express Gaussian Elimination formally as 

(**U, c) \= GN**\-1PN-**IGN**\-2PN-*2* G1P1 (**A, b**).   
···   
**(A,b)** 

\=   
\[16\] 

With these definitions we can touch upon another very important class of direct methods, that of LU decompostion. The idea is to decompose the original matrix **A** into a product of an upper triangular matrix U and a lower triangular matrix L, i.e. **A LU**. Such decomposition is important because, once we have it, we can solve the system of equations for may right hand sides. This is easy seen. The original system was **Ax** \= **b.** With the LU decomposition this becomes **LUX** \= b. First solve **Ly** \= **b** using forward substitution (the same trick as with backward substitution, but now starting with y1 and working in the forward direction). Next, solve **Ux \= y** through backward substitution.   
\=   
\=   
*Numerical Methods*   
17 

In the special case of Gaussian Elimination *without* the need for row exchanges, i.e. when P; \= I, we can derive a **LU** decomposition using the multipliers li. In this special case, from equation \[16\] we can immediately derive that 

**U**\=GN-**IGN**\-2G, **A ALU**, where L=GGG. Furthermore, it is easily verified that 

1   
0 

**L**\=   
:   
1 0 

INI   
*IN*,N-I 1 

3   
1 6   
X1   
23 

As an example, consider   
1   
13 X2   
12   
With Gaussian Eliminiation one can 

1 1 1 *X3*   
6 

easily verify that   
X2 

X3   
0-0   
1   
0   
0   
3 1 6 

**L**\= 1/3 1   
9 

3   
0 and **U**\=0 2/3 1 

1/3 1 1   
0 0 \-2   
As an 

exercise perform the calculation yourself and check the results. 

...   
1   
In general **LU** \= **PA** with **P** \= PN-1 PN-2 P1 and L some permutation of elements from the **G**; matrices. Many algorithms exist to find LU factorizations (see e.g. \[2\]). Finally note that also many optimized algorithms exist for special matrices, e.g for banded matrices. For more information on LU decomposition read section **Error\! Reference source not found.** and **Error\! Reference source not found.** of the appendix. 

Parallel Libraries for Direct Methods 

The scientific computing community has invested a lot of time in creating very powerful and highly efficient (parallel) public domain libraries for numerical algorithms. Many of them can be found on the famous *netlib* web site,.3 which you are advised to visit and use the available software whenever needed. 

In most linear algebra algorithms one typically encounters a small set of basic routines. These basic routines are formalized into the BLAS set (Basis Linear Algebra Subroutines, for more in-depth information, see the appendix, **Error\! Reference source not found..** The BLAS set is divided into three levels, 

BLAS 1: O(N) operations on O(N) data items, typically a vector update, 

e.g. : **y** \= ax+**y**; 

BLAS 2: O(N2) operations on O(N2) data items, typically matrix vector products, e.g. **: y** \= **aA x** *\+b* **y**; 

BLAS 3: O(N3) operations on O(N2) data items, typically matrix-matrix products, e.g. **CaA B** *\+b* **C.** 

Highly optimized BLAS libraries are available on e.g. netlib. Performance increases in going from BLAS-1 to BLAS-3 due to better ratio between amount of calculations and needed memory references. BLAS3 was specifically developed for use on parallel computers. Gaussian elimination or LU factorization can be implemented using BLAS-1, BLAS-2, or BLAS-3. Modern linear algebra libraries use BLAS-3, where the algorithms are formulated as block algorithms. An example of such a block LU factorization can be found in the appendix, section **Error\! Reference source not found..** 

Using BLAS routines a number of well-known linear algebra libraries have been developed and put into the public domain. As examples we mention: 

3 www.netlib.org   
*Numerical Methods*   
18 

**D.3.**   
EISPACK, Fortran routines for calculation of eigenvectors and eigenvalues for dense and 

banded matrices, based on BLAS-1; 

LINPACK, Fortran routines to solve linear equations for dense and banded matrices, 

column oriented BLAS-1 approach; 

LAPACK, successor of EISPACK and LINPACK, Fortran routines for eigenvalues and solving linear equation, for dense and banded matrices, exploit BLAS-3 (to improve speed), runs efficiently on vector machines and distributed memory parallel computers; 

ScaLAPACK, extend LAPACK to run efficiently on distributed memory parallel 

computers, fundamental building blocks are distributed memory versions of BLAS-2 and BLAS-3 routines and a set of Basic Linear Algebra Communication Subprograms (BLACS). 

We end this section with noting that a detailed discussion of parallel algorithms for e.g. LU factorization is beyond the scope of this lecture. The main features of such parallel algorithms are a block oriented BLAS-3 approach in combination with special (scattered) decompositions of the matrices. For more information we refer to e.g. Reference \[4\], chapter 9.5 or Reference \[5\]. 

**Iterative Methods** 

The General Idea 

\= 

...   
\-   
As was already explained before, in iterative methods one seeks an approximate solution via some iterative procedure. The idea is to consider an iteration of the form x(\+1) *P(*x("), n 0, 1, 2, The superscript n is the iteration number. Given a first guess solution x and the iteration function one iterates the solution forward until some convergence criterion is met. Typically one demands that ||x(\+1) x(")|| \< *ε,* where *ε* is some small number. 

For our linear system **Ax** note that   
\=   
**b** an iterative procedure can be constructed as follows. First 

**Ax** \= **b → Bx** \+ **(A \- B)x** \= **b** 

then put 

***By***(n*\+1\)* 

**or**   
\+ **(A \- B)x**(") \= b 

x(n\+1) \=   
\= (I**\-B1A)x**(") \- **B'b** 

\-   
Many possibilities exist for the choice of the matrix **B**, thus leading to many different types of iterative methods (see e.g. \[2\], chapter 8\). In this section we will only consider three more or less related iterative methods, the *Jacobi* iteration, the *Gauss-Seidel* iteration, and *Successive Over Relaxation*. Reference \[2\] and the appendix, sections **Error\! Reference source not found. Error\! Reference source not found.,** provide many more examples of iterative methods. Here we especially mention the family of iterative methods called *Krylov space* methods, also knows as *Conjugate Gradient* methods. They are very powerful, in the sense that they require not so many iterations and that they can handle bad conditioned systems as well. In many cases they are *the* method of choice as a solver. Like the iterative methods that will be introduced here they are relatively easy to parallelize. Finally, many efficient parallel implementations of iterative methods are available in the public domain and a highly readable book that provides a cookbook style of using iterative methods (Reference \[6\]) can be obtained in electronic from netlib.   
*Numerical Methods*   
19 

The Jacobi Iterative Method 

Let us now return our attention to the specific case of the finite difference scheme for the two\-dimensional Laplace equation (Eq. \[14\]). This equation immediately suggests an iterative scheme: 

*(n*\+1\) ci,m   
\=   
*(**n***) Citi,m 4   
*(n*) \+cm   
*(n*) \+cm\+1 \+ *Ci\=1*,m \+ *Ci,m*\+1 \+c   
*(**n***) *Ci,m***\-1**   
\[17\] 

As before, *n* is the iteration index. This iterative scheme is known as the Jacobi iteration. Note the relation of Equation \[17\] with the finite difference scheme for the time dependent diffusion equation, Eq. \[11\]. If we take the by the CFL condition allowed maximum time step, i.e. *D\&x2/\&t* 1/4 and use that in Eq. \[11\] we reproduce the Jacobi iteration of Eq. \[17\]. In other words, we expect that the Jacobi iteration will also suffer from a relative large amount of iterations needed before convergence.   
\= 

Algorithm 2 provides the code for the inner loops of a Jacobi iteration. We again assume the square domain with periodic boundaries in the x\-direction and fixed boundaries in the y-direction. Furthermore, a specific stopping criterion is implemented. We demand that 

max *ij*   
*(n*\+1\)   
' \- c}") | \< *ɛ* . 

The difference in concentration between two iterations on all grid points should be smaller that some small number *ε*. This is a rather severe stopping condition. Others can also be used, e.g. calculating the mean difference and demanding that this should be smaller that some small number. Here we will not pay any further attention to the stopping criterion, but you should realize this is an important issue to be considered in any new application of an iterative method. 

/\* Jacobi update, square domain, periodic in x, fixed \* /\* upper and lower boundaries 

**do** {   
\*   
/ 

δ   
\=   
0 

**for** i\=0 to max { 

**for** j\=0 to max { 

*(n\+*1)   
\=   
1.0 

*(n\+*1)   
\=   
0.0   
**if**(*cij* is a source) *Cij* **else if** (*Cij* is a sink) *Cij* **else** { 

/\* periodic boundaries \*/ 

west 

east   
\= 

\=   
(n)   
(*n*) 

(*n*)   
(n)   
: Ci\+l, j   
(i\==0\) ? Cmax\-1,j : *Ci\-*1,j (i==max) ? C1,j 

/\* fixed boundaries \* / 

Cij   
south 

north \= 

*(n*\+1) \=   
\=   
(j==0) ? c0 :   
*(n)*   
Ci,j\-1 

(j==max) ? cL 

0.25 \*   
(west   
(*n)*   
:   
Ci, j\+1 

} 

/\* stopping criterion \* /   
east \+ south \+ north) 

**if** (*Cij*   
(*n*\+1)   
—   
*C(*) \> tolerance) *8* \=   
*(n*\+1)   
(*n*)   
\=   
Cij   
*Cij*   
ן (n 

} 

} 

while (8\> tolerance) 

*Algorithm 2: The sequential Jacobi iteration* 

As an example again we take the two\-dimensional square domain, 0 ≤ x, *y* ≤ 1 with periodic boundary conditions in x-direction, *c*(x,y) \= c(x+1,y), and fixed values for the upper and lower boundaries, c(x,y\=0\) \= 0 and *c(x*,y\=1\) \= 1\. Because of symmetry the   
*Numerical Methods*   
20   
20 

solution will not depend on the *x*\-coordinate, and that the exact solution is a simple linear concentration profile: *c(x*,y) \= *y* (derive this result yourself\!). The availability of the exact solution allows us to measure the error in an iterative method as a function of the applied stopping criterion. We will show some results after having introduced the other iterative methods, allowing for a comparison. Here we first concentrate on the number of iterations needed for convergence. 

Set the small number *&* in the stop condition to *ɛ*\= 107, where *p* is a positive integer. Next we measure the number of iterations that is needed for convergence, and we do this for two grid sizes, *N* \= 40 and *N* \= 80 (where *N* \= *L* \+ 1). As first guess for the solution we just take all the concentration equal to zero. The results, shown in Figure 12, seem to suggest a linear dependence between the number iteration and *p*, i.e. a linear dependence on \-log(*)*. Furthermore, the results suggest an *N2* dependence. These results can be obtained from mathematical analysis of the algorithm, see e.g. \[2, 3\]. This *N2* dependence was also found for the time dependent case (see discussion in section C.1). Finally note that the number of iterations can easily become much larger than the number of grid points *(*N) and in that case direct methods to solve the equations are preferred. In conclusion, the Jacobi iteration works, but only as an example and certainly not to be used in real applications. However, the Jacobi iteration is a stepping stone towards a very efficient iterative procedure.   
*iterations*   
3.0E+04 

2.5E+04 

2.0E+04 

1.5E+04 

1.0E+04 

5.0E+03 

0.0E+00 

0   
2 

*р*   
10   
N \=40 

N \= 80 

*Figure 12: The number of iterations for the Jacobi iteration as a function of the stopping condition and as a function of the number of grid points along each dimension*. 

At first sight a parallel Jacobi iteration seems very straightforward. The computation is again based on a local five-point stencil, as in the time dependent case. Therefore, we can apply a domain decomposition, and resolve all dependencies by first exchanging all boundary points, followed by the update using the Jacobi iteration. However, the main new issue is the stopping condition. In the time dependent case we know before run time how many time steps need to be taken. Therefore each processor knows in advance how many iterations are needed. In the case of the Jacobi iteration (and all other iterative methods) we can only decide when to stop during the iterations. We must therefore implement a parallel stopping criterion. This posses a problem, because the stopping condition is based on some global measure (e.g. a maximum or mean error over the complete grid). That means that a global communication (using e.g. MPI-Reduce) is needed that may induce a large communication overhead. In practical implementations one must seriously think about this. Maybe it pays off to calculate the stopping condition once every q iterations (with q some small positive number) instead of after each iteration. This depends on many details, and you are challenged to think about this yourself (and apply your ideas in the lab course that is associated to this lecture). With all this in mind the pseudo code for the parallel Jacobi iteration is given in Algorithm 3\. It is clear that it closely resembles the time-dependent pseudo code.   
*Numerical Methods* 

**main** () /\* pseudo code for parallel Jacobi iteration 

{ 

*decompose lattice*;   
\*   
21 

}   
*initialize lattice sites;* 

set *boundary conditions;* 

**do** { 

}   
*exchange boundary strips* with *neighboring processors;* **for** *all grid points in this processor* { 

}   
update according to Jacobi iteration; 

calculate local *&* parameter; /\* stopping criterion 

obtain the global maximum *♂* of all local *&* values 

**while** (*8*\> tolerance) 

*print results to file*; 

*Algorithm 3: Pseudo code for the parallel Jacobi iteration*. 

The Gauss-Seidel Iterative Method   
\*   
/ 

The Jacobi iteration is not very efficient, and here we will introduce a first step to improve the method. The *Gauss*\-*Seidel* iteration is obtained by applying a simple idea. In the Jacobi iteration we always use results from the previous iteration to update a point, *even* when we already have new results available. The idea of Gauss-Seidel iteration is to apply new results as soon as they become available. In order to write down a formula for the Gauss-Seidel iteration we must specify the order in which we update the grid points. Assuming a row-wise update procedure (i.e. we increment / while keeping *m* fixed) we find for the Gauss-Seidel iteration 

*(n*\+1\) ci,m   
\=   
H   
*(n) Ci\+*1,m   
\~ (n \+ 1) \+ C { m2 \+ 1   
*(n*) \+ C1-1,m   
(n\+1) \+ C1.m\-1   
4   
\[18\] 

One immediate advantage of the Gauss-Seidel iteration lies in the memory usage. In the Jacobi iteration you would need two arrays, one to store the old results, and another to store the new results. In Gauss-Seidel you immediately use the new results as soon as they are available. So, we only need one array to store the results. Especially for large grids this can amount to enormous savings in memory\! We say that the Gauss\-Seidel iteration can be computed *in place*. 

Is Gauss-Seidel iteration also faster than Jacobi iteration. According to theory it turns out that a Gauss\-Seidel iteration requires a factor of two iterations less then Jacobi (see \[3\], section 17.5). This is also suggested by a numerical experiment. We have taken the same case as in the previous section, and for *N* \= 40 we have measured the number of iterations needed for Gauss-Seidel and compared to Jacobi. The results are shown in Figure 13\. The reduction of the number of iterations with a constant number is indeed observed, and this constant number is very close to the factor of two as predicted by the theory. This means that Gauss-Seidel iteration is still not a very efficient iterative procedure (as compared to direct methods). However, Gauss-Seidel is also only a stepping stone towards the Succesive Over Relaxation method (see next section) which *is* a very efficient iterative method. 

The Gauss\-Seidel iteration posses a next challenge to parallel computation. At first sight we must conclude that the parallelism available in Jacobi iteration is now completely destroyed by the Gauss-Seidel iteration. Gauss\-Seidel iteration seems inherently sequential. Well, it is in the way we introduced it, with the row-wise ordering of the computations. However, this row\-wise ordering was just a convenient choice. It turns out that if we take another ordering of the computations we can restore parallelism in the   
*Numerical Methods*   
**22**   
22 

Gauss-Seidel iteration. This is an interesting case where reordering of computations provides parallelism. Keep this in mind, as it may help you in the future in finding parallelism in algorithms\!   
*iterations*   
8.0E+03 

7.0E+03 

6.0E+03 

5.0E+03 

4.0E+03 

3.0E+03 

2.0E+03 

1.0E+03 

0.0E+00 

0   
**2**   
6   
8   
10 

*p*   
Jacobi 

Gauss-Seidel 

*Figure 13: The number of iterations for the Jacobi and Gauss*\-*Seidel iteration as a function of the stopping condition for* N \= *40*. 

The idea of the reordering of the computations is as follows. First, color the computational grid as a checkerboard, with red and black grid points. Next, given the fact that the stencil in the update procedure only extends to the nearest neighbors, it turns out that all red points are independent from each other (they only depend on black points) and vice-versa. So, instead of the row-wise ordering we could do a red-black ordering, were we first update all red points, and next the black points (see Figure 14). We also call this Gauss\-Seidel iteration, because although the order in which grid points are updated is now different, we also do the computation in place, and use new results as soon as they become available. 

*Figure 14: Row wise ordering (left) versus red-black ordering (right)*. 

This new red-black ordering restores parallelism. We can now first update all red points in parallel, followed by a parallel update of all black point. The pseudo code for parallel Gauss-Seidel with red-black ordering is show in Algorithm 4\. The computation is now split into two parts, and before each part a communication with neighboring processors is needed. On first sight this would suggest that the parallel Gauss-Seidel iteration requires twice as many communication time in the exchange part. This however is not true, because it is not necessary to exchange the complete set of boundary points, but only half. This is because for updating red points we only need the black point, so also only the black boundary points need to be exchanged. This means that, in comparison with parallel Jacobi, we only double the setup times required for the exchange operations, but keep the sending times constant. Finally note that the same parallel stop condition as before can be applied.   
*Numerical Methods* 

/\* only the inner loop of the parallel Gauss\-Seidel method with /\* Red Black ordering \*/ 

**do** { 

*exchange boundary strips with neighboring processors;* **for** *all* red *grid points in this processor* { 

}   
update according to Gauss\-Seidel iteration; 

*exchange boundary strips with neighboring processors;* 

**for** *all* black *grid points in this processor* {   
\* 

update according to Gauss-Seidel iteration; 

} 

obtain the global maximum *8* of all local *8* values 

}   
23 

**while** (*8*\> tolerance) 

*Algorithm 4: The pseudo code for parallel Gauss-Seidel iteration with red-black ordering*. 

Another final remark is that instead of applying red-black ordering to each grid point in the domain, we can also create a coarse\-grained red\-black ordering, as drawn in Figure 15\. Here the red-black ordering is done on the level of the decomposed grid. The procedure could now be to first update the red domain followed by an update of the black doamin. Within each domain updating can now be done with the row-wise ordering scheme (why?). 

Ро   
1   
P1 P2 

*Figure 15: A coarse\-grained red-black ordering*. 

Successive Over Relaxation 

The ultimate step in the series from Jacobi and Gauss\-Seidel is to apply a final and as will become clear very efficient idea. In the Gauss-Seidel iteration the new iteration results is completely determined by its four neighbors. In the final method we apply a correction to that, by mixing the Gauss-Seidel result with the current value, i.e. 

c{(n \+ 1) \= @ \[c (n) \+ *c*{(n \+ 1\) \+ c{(m)\+ \+ *c*{\+} \] \+ (1 − *@) c (* ).   
(*n*\+1)   
(*n*\+1\)   
Ci,m   
\] \- *0\)* (   
4   
(n\+1\) Ci-1,m   
\-C12m\+1   
·   
\[19\] 

The parameter *@* determines the strength of the mixing. One can prove (see e.g. \[2\]) that for 0 \< *\<* 2 the method is convergent. For *@* \= 1 we recover the Gauss-Seidel iteration, for 0 \< \< 1 the method is called Successive Under Relaxation. For 1 \< *@*\< 2 we speak of *Successive Over Relaxation,* or SOR. 

It turns out that for finite difference schemes SOR is very advantageous and gives much faster convergence then Gauss-Seidel. The number of needed iterations will however strongly depend on the value of *@*. One needs to do some experiments in order to determine its optimum value. In our examples we have take *@* always close to 1.9. The enormous improvement of SOR as compared to Gauss-Seidel and Jacobi is illustrated by   
*Numerical Methods* 

again measuring the number of iterations for the example of the square domain with *N* \=40\. The results are shown in Figure 16, and show the dramatic improvement of SOR. Another important result is the final accuracy that is reached. Remember that we know the exact solution in our example. So, we can compare the simulated results with the exact solution and from that calculate the error. We have done that and the results are shown in Figure 17\. Note the logarithmic scale on error-axis in Figure 17\. Because we observe a linear relationships between the logarithm of the error and *p* we can conclude that we find a linear relationship between the stopping condition *ɛ* and the mean error in the simulations. Furthermore we observe that the error in the SOR is much smaller than in Jacobi and Gauss-Seidel. This further improves the efficiency of SOR. Suppose you would like to get mean errors of 0.1%. In that case in SOR we only need to take *p* 4, whereas Jacobi or Gauss-Seidel require *p* \= 6\. SOR is a practically useful iterative method and as such is applied regularly. Finally note that parallel SOR is identical to parallel Gauss\-Seidel. Only the update rules have been changed.   
\=   
*iterations*   
8.0E+03 

7.0E+03 

6.0E+03 

5.0E+03 

4.0E+03 

3.0E+03 

2.0E+03 

1.0E\+03 

0.0E+00 

0   
2   
8   
10 

\- Jacobi   
\-Gauss-Seidel   
SOR 

*Figure 16: The number of iterations for the Jacobi, Gauss\-Seidel, and SOR iteration as a function of the stopping condition for* N \= *40\.*   
*error*   
1.0E+00 

1.0E-01 

1.0E-02 

1.0E-03 

1.0E-04 

1.0E-05 

1.0E-06 

1.0E-07 

0   
2   
8   
10 

*p* 

Jacobi   
Gauss-Seidel SOR 

*Figure 17: The mean error for the Jacobi, Gauss-Seidel, and SOR iteration as a function of the stopping condition for* N \= *40*. 

Iterative Methods in Matrix Notation 

So far we only considered the iterative methods in the special case of the finite difference discretization of the two-dimensional Laplace equation. Now we return to the general idea   
24   
*Numerical Methods*   
25   
25 

of constructing iterative methods. Remember that in general an iterative method can be constructed from **x'**   
(*n*\+1)   
(I \- **B1A)x("**) **\- B ́1b.** Now with   
\= 

0   
a11 

**D=**   
**A**\=D+**E+F**, 

0   
a NN 

0   
0 a12   
*a1N*   
0 

**E** \=   
a21 

:   
**F**\= 

a *N*\-1*,N* 0   
0   
*a* N1   
a *N*,*N*\-1 

we can define in a general way the three iterative methods introduced so far. 

Jacobi iteration 

**B\=**D \= *ax*(n+1)\=\-Σ *a ¡¡ x*("*) \+b*; 

i\+*j* 

Gauss-Seidel iteration 

**B**\=**D\+E** ⇒ α *¡¡* x ;   
*(n*\+1)   
\== 

SOR   
(n \+ 1) \= \[ *a1* x (n*) \+ b* ; 

*i\>j*   
\-Σ a *1j* x} 

i\<*j* 

\=   
**B** \-D**\+E** ⇒aji *x* j   
*(n*\+1)   
\=   
***(n***)   
\- Σ a ¡¡x (\+1) \- Σ a ¡¡x\!") \+b;) \+(1−*@*)*a;¡¡x{"*)   
\-Σ 

*i\<j*   
*i\>j* 

**D.4.**   
**An Application: Diffusion Limited Aggregation** Diffusion Limited Aggregation (DLA) is a model for non-equilibrium growth, where growth is determined by diffusing particles. It can model e.g. a Bacillus subtilis bacteria colony in a petri dish. The idea is that the colony feeds on nutrients in the immediate environment, that the probability of growth is determined by the concentration of nutrients and finally that the concentration of nutrients in its turn is determined by diffusion. The basic algorithm is shown in Algorithm 5. 

1\. Solve Laplace equation to get distribution of nutrients, assume that the object is a sink (i.e. *c* \= 0 on the object) 2\. Let the object grow 

3\. Go back to (1) 

*Algorithm 5*: *The DLA algorithm*. 

The first step in Algorithm 5 is done by a parallel SOR iteration. Step 2, growing of the object, requires three steps; 

1.   
determine growth candidates; 

2\. determine growth probabilities; 

3\.   
grow. 

A growth candidate is basically a lattice site that is not part of the object, but whose north \-, or east   
or south \-, or west neighbor is part of the object. In Figure 18 a possible configuration of the object is shown; the black circles form the current object, while the white circles are the growth candidates.   
*Numerical Methods* 

*Figure 18: The object and possible growth sites*. 

The probability for growth at each of the growth candidates is calculated by 

°   
*Pg ((i*, *j)* \= • → *(*i, j)€ •)=   
(*ci,j)n* Σ(ci,;*)"*   
\[20\] 

*(i, j*)Є° 

\=   
The parameter *n* determines the shape of the object. For *n* 1 we get the normal DLA cluster, i.e. a fractal object. For *n* \< 1 the object becomes more compact (with ŋ \= 0 resulting in the Eden cluster), and for *ŋ* 1 the cluster becomes more open (and finally resembles say a lightning flash). 

Modeling the growth is now a simple procedure. For each growth candidate a random number between zero and one is drawn and if the random number is smaller than the growth probability, this specific site is successful and is added to the object. In this way, on average just one single site is added to the object. 

The results of a simulation are shown in Figure 19\. The simulations were performed on a 2562 lattice. SOR was used to solve the Laplace equation. As a starting point for the SOR iteration we used the previously calculated concentration field. This reduced the number of iterations by a large factor. For standard DLA growth *(n* \= 1.0) we obtain the typical fractal pattern. For Eden growth a very compact growth form is obtained, and for *n* \= 2, a sharp lighting type of pattern is obtained. 

100 

50   
250 

200   
250 

200 

**150**   
**150**   
150 

50   
100   
150   
100   
100   
250 

100   
150   
250 

*Figure 19: Results of DLA growth on a 2562 lattice*. *The left figure is for n* \= *0*, *the middle for n 1.0 and the right for n \= 2.0*.   
\= 

Introduction of parallelism in DLA posses an interesting new problem (that we will not solve here, but in Chapter 4\). It is easiest is to keep using the same decomposition as for the iterative solvers (e.g. the strip-wise decomposition of the grid). The growth step is local, and therefore can be computed completely in parallel. Calculation of the growth probabilities requires a global communication (for the normalization, i.e. the denominator of the equation for the probability). However, due to the growing object the load balancing gets worse and worse during the simulation. Our computational domain is *dynamically* changing during the simulation. This calls for a completely different view on decomposition, and that will be discussed in detail in Chapter 4.   
26   
26 