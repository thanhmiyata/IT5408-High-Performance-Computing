Proceedings of the Czech–Japanese Seminar in Applied Mathematics 2006 

Czech Technical University in Prague, September 14-17, 2006   
pp. 25–36 

PARALLEL ALGORITHM FOR NUMERICAL SOLUTION OF THE SHALLOW WATER EQUATION 

STANISLAV BRAND1 

Abstract. Parallel programming is the leading technique for accelerating numerical algorithm. This article deals with parallelization of the program for solving the shallow water equation in 2D, representing an example of a conservation law. The equation is solved by finite difference and finite volume methods within the Lax-Friedrichs, Lax-Wendroff and MacCormack scheme. The algorithm is parallelized using OpenMP and MPI. The efficiency measurement is made mainly on systems with shared memory and thanks to the sufficient number of processors there are some results of the mixed mode programming which used both OpenMP and MPI. 

Key words. Conservation laws, parallel programming, OpenMP, MPI, mixed mode programming. 

1\. Introduction. This article summarizes the results of numerical solution of the shal low water equation. The first section shows the derivation of the one-dimensional shallow water equation and describes the two-dimensional shallow water system. The second section contains numerical schemes used, in the particular finite difference and finite volume vari ants of the well-known Lax-Friedrichs, Lax-Wendroff and MacCormack schemes. The third section describes parallelization, that was made for improving program performance. The fourth and fifth section contains several results of efficiency measurement and of the problem solved. 

2\. Solved problem. The Shallow water equation is a special case of the Euler equations describing wave motion in shallow water. The two-dimensional shallow water equation is often used for numerical simulation of natural river flows. For illustration, we derive the one dimensional shallow water equation which is also often used. Consider a fluid in a channel and assume that the vertical velocity in the fluid is negligible and the horizontal velocity *v*(*x, t*) is constant through any vertical cross section. For example, small aptitude wave in the fluid that is shallow relative to the wave length. We assume an incompressible fluid, where the density *ρ* is constant. Our variable is the height of water *h*(*x, t*). The total mass in *x*1*, x*2   
at time *t* can be written as *mx*1*,x*2 \=R *x*2   
*x*1*ρh*(*x, t*)d*x*. Vertical integration from 0 to *h*(*x, t*)   
of the momentum *ρv*(*x, t*) gives us the mass flux to be *ρv*(*x, t*)*h*(*x, t*). By dropping out the constant *ρ* we get the conservation of mass equation in the form 

*ht* \+ (*vh*)*x* \= 0*.* (2.1) 

Now we get the conservation of momentum from the Euler equation 

(*ρhv*)*t* \+ (*ρhv*2 \+ *p*)*x* \= 0*.* (2.2) 

In our case, the pressure *p* is determined from a hydrostatic law saying the pressure at depth *y* is *ρgy*, where *g* is the gravitation constant. Integrating this relation vertically from *y* \= 0 to *y* \= *h*(*x, t*) we get the proper pressure term *p* \=12*ρgh*2. Using this in (2.2), and pointing out *ρ* gives 

1Department of Mathematics, Faculty of Nuclear Sciences and Physical Engineering, Czech Technical Univer sity, Prague. 

25  
26 S. Brand 

z 

b(x,y) 

y 

x 

h(x,y) 

FIG. 2.1. *Meaning of variables for the shallow water equation.* 

(*hv*)*t* \+ (*hv*2 \+12*gh*2)*x* \= 0*.* (2.3) 

One-dimensional shallow water system consists of the continuity equation (2.1) and the momentum equation (2.3). We use the two-dimensional shallow water equation whose deriva tion can be found in \[3\]. 

It can be written as 

*Ut* \+ *F*(*U*)*x* \+ *G*(*U*)*y* \= *S,* (2.4) 

where 

*U* \=   
   
*hhu* 

 

 \= 

   
*u*1 *u*2   
 

 *,* 

   
*hu* 

*hv*   
 

   
*u*3 

*u*2 

 

*F*(*U*) \= *G*(*U*) \=   
 

    
*hu*2 \+12*gh*2 *huv* 

*hv* 

*huv* 

*hv*2 \+12*gh*2    
 \=    
 \= 

 

    
*u*22*/u*1 \+12*gu*21 *u*2*u*3 

*u*1 

*u*3   
*u*2*u*3   
*u*1   
*u*23*/u*1 \+12*gu*21 

 *,*    
 *,* 

*S* \=   
0   
*−ghbx*   
*−ghby*   
 

 *.* 

Here *h* \= *h*(*x, y*) denotes the depth of water over the bottom topography *b* \= *b*(*x, y*), *u*,*v* are velocities in *x* and *y* directions, and *g* stands for the gravitational constant (see Figure 2.1). 

3\. Numerical solution. We used three finite difference methods for conservation laws, the Lax-Friedrichs(LF), the Lax-Wendroff(LW), and the MacCormack(MC) schemes. For details, we refer the reader to \[4\],\[5\]. All of these methods have some limitations. The LF scheme is first-order-accurate and diffusive. The LW and MC schemes are second-order accurate, but oscillatory near shocks.  
Parallel algorithm for numerical solution of problems in fluid dynamics 27 

We therefore also tested the composite scheme as well which combines these methods in order to eliminate their problems. For more information about the composite schemes, see \[6\],\[7\]. All computations were made with the composite scheme that was a composition of three LW steps and one diffusive LF step needed for removing spurious oscillations. 

The schemes in the following sections are based on the finite volume and finite difference discretization. 

3.1. Lax-Friedrichs scheme. Finite volume method yields: 

*Fni,k*∆*yk − Gni,k*∆*xk,* (3.1) 

where   
*i* \=1*N*X*N*   
*Un*\+1 

*k*\=1   
*Unk −*∆*t µ*(*Di*)   
X*N k*\=1 

*Fni,k* \=12\[*F*(*Uni*) \+ *F*(*Unk*)\]   
*Gni,k* \=12\[*G*(*Uni*) \+ *G*(*Unk*)\]*.*(3.2) 

Finite difference method yields: 

*i,j* \=14(*Uni,j−*1 \+ *Uni,j*\+1 \+ *Uni−*1*,j* \+ *Uni*\+1*,j* )   
*Un*\+1   
*−k*2*h*£*F*(*Uni*\+1*,j* ) *− F*(*Uni−*1*,j* )¤*−k*2*l*£*G*(*Uni,j*\+1) *− G*(*Uni,j−*1)¤*.*(3.3) 3.2. Lax-Wendroff scheme. Finite volume method yields: 

*i* \=1*N*X*N*   
*Un*\+ 12 

*k*\=1 

*Unk −*∆*t* 2*µ*(*Di*)   
X*N k*\=1   
*Fni,k*∆*yk − Gni,k*∆*xk,* 

(3.4)   
*i* \= *Uni −*∆*t*   
*Un*\+1 

*µ*(*Di*) 

where   
X*N k*\=1 

*Fn*\+ 12   
*i,k* ∆*yk − Gn*\+ 12   
*i,k* ∆*xk,* 

*i,k* \=12\[*F*(*Un*\+ 12   
*Fni,k* \=12\[*F*(*Uni*) \+ *F*(*Unk*)\]*, F n*\+ 12 

*i*) \+ *F*(*Un*\+ 12 *k*)\]*,*   
(3.5)   
*i,k* \=12\[*G*(*Un*\+ 12   
*Gni,k* \=12\[*G*(*Uni*) \+ *G*(*Unk*)\]*, Gn*\+ 12 

Finite-difference method yields: 

*i,j* \=14(*Uni,j−*1 \+ *Uni,j*\+1 \+ *Uni−*1*,j* \+ *Uni*\+1*,j* ) *Un*\+ 12 

*i*) \+ *G*(*Un*\+ 12 *k*)*.* 

*Un*\+1   
*−k*4*h*\[(*F*(*Uni*\+1*,j* ) *− F*(*Uni−*1*,j* )\] *−k*4*l*\[(*G*(*Uni,j*\+1) *− G*(*Uni,j−*1)\] 

(3.6)  
*i,j* \= *Uni,j*   
*−k*2*h*\[(*F*(*Un*\+ 12 

*i−*1*,j* )\] *−k*2*l*\[(*G*(*Un*\+ 12   
*i*\+1*,j* ) *− F*(*Un*\+ 12   
*i,j*\+1) *− G*(*Un*\+ 12 *i,j−*1)\]*.*   
28 S. Brand 3.3. MacCormack schemes. Finite volume method yields: 

*i* \= *Uni −*∆*t*   
X*N*   
*Un*\+ 12 

" 

*µ*(*Di*)   
*Fni,k*∆*yk − Gni,k*∆*xk, k*\=1 

\#   
(3.7)   
*i* \=12   
*Un*\+1 

where 

*i −*∆*t*   
*Uni* \+ *Un*\+ 12 

*µ*(*Di*)   
X*N k*\=1   
*Fn*\+ 12   
*i,k* ∆*yk − Gn*\+ 12 

*i,k* ∆*xk*   
*,* 

*Fni,*1 \= *F*(*Uni*\+1*,j* )*, F ni,*2 \= *F*(*Uni,j*\+1)*, F ni,*3 \= *Fni,*4 \= *F*(*Uni,j* ) *Fn*\+ 12   
*i,*3 \= *F*(*Un*\+ 12   
*i−*1*,j* )*, F n*\+ 12   
*i,*4 \= *F*(*Un*\+ 12 

*i,j−*1)*, F n*\+ 12   
*i,*1 \= *Fn*\+ 12   
*i,*2 \= *F*(*Un*\+ 12   
*i,j* ) 

*Gni,*1 \= *G*(*Uni*\+1*,j* )*, Gni,*2 \= *G*(*Uni,j*\+1)*, Gni,*3 \= *Gni,*4 \= *G*(*Uni,j* ) *Gn*\+ 12   
*i,*3 \= *G*(*Un*\+ 12   
*i−*1*,j* )*, Gn*\+ 12   
*i,*4 \= *G*(*Un*\+ 12 

*i,j−*1)*, Gn*\+ 12   
*i,*1 \= *Gn*\+ 12   
*i,*2 \= *G*(*Un*\+ 12   
*i,j* )*.* 

Finite difference method yields: 

*Un*\+ 12   
*i,j* \= *Uni,j*   
*−kh*\[(*F*(*Uni*\+1*,j* ) *− f*(*Uni,j* )\] *−kh*\[(*G*(*Uni,j*\+1) *− G*(*Uni,j* )\]   
*i,j* \=12(*Uni,j* \+ *Un*\+ 12 

(3.8) (3.9)   
*Un*\+1 

*i,j* ) 

*−k*2*h*\[(*F*(*Un*\+ 12 

*i−*1*,j* )\] *−k*2*h*\[(*G*(*Un*\+ 12   
*i,j* ) *− F*(*Un*\+ 12   
*i,j* ) *− G*(*Un*\+ 12 *i,j−*1)\]*.* 

4\. Stability. Each numerical scheme must satisfy the necessary stability condition in the form 

∆*t ≤*1 ∆*~~x~~* \+*~~ρ~~B*   
*,* (4.1)   
*~~ρ~~A*   
∆*y* 

where *ρA* and *ρA* are spectral radii of jacobian matrix *A* and *B*,    
 

*A* \=*∂F*(*U*) *∂U* \=   
0 1 0 *gh − u*2 2*u* 0 

 *,* (4.2) 

   
*−uv v u* 

 

*B* \=*∂G*(*U*) *∂U* \=   
0 0 1   
*−uv v u gh − v*2 0 2*v* 

 *.* (4.3) 

Unfortunately, this condition is not sufficient for nonlinear partial differential equations. To demonstrate stability and consistency, we have to use numerical results computed on a refined grid. We linearly interpolate the solution on the finest grid and compare it with the remaining solutions (see Tables 4.1 \- 4.4). 

5\. Numerical results. This section contains results for the shallow water equation with the flat bottom. Results are represented by graphs displaying depth of water *h* at each point of the space domain (see Figures 5.1 \- 5.3). The initial condition for the depth of water  
Parallel algorithm for numerical solution of problems in fluid dynamics 29 

| Mesh  *h*  | *L∞*(0*, T*;*L*2) error of *u*  | *L∞*(0*, T*;*L∞*)error of *u* |
| :---- | :---- | :---- |
| 0.0800000 0.0400000 0.0200000 0.0100000 0.0050000  | 0.0154219  0.0078344  0.0036292  0.0015629  0.0005506  | 0.0149480 0.0071600  0.0033470  0.0016555  0.0006410 |

T~~ABLE~~ 4.1  
*Table of convergence errors for the MacCormack scheme.* 

| Mesh  *h*  | EOC *u  L*2  | EOC *u L∞* |
| :---- | :---- | :---- |
| 0.0800000 0.0400000 0.0200000 0.0100000 0.0050000  | 0.0000000 0.9584786 1.3806293 0.9930047 1.0094179  | 0.00000001.2418251 1.0897255 1.1842024 1.2039963 |

T~~ABLE~~ 4.2  
*Table of EOC coefficients for the MacCormack scheme.* 

| Mesh  *h*  | *L∞*(0*, T*;*L*2) error of *u*  | *L∞*(0*, T*;*L∞*)error of *u* |
| :---- | :---- | :---- |
| 0.0800000 0.0400000 0.0200000 0.0100000 0.0050000  | 0.0154219  0.0078344  0.0036292  0.0015629  0.0005506  | 0.0149480 0.0071600  0.0033470  0.0016555  0.0006410 |

T~~ABLE~~ 4.3  
*Table of convergence errors for the composite scheme.* 

| Mesh  *h*  | EOC *u  L*2  | EOC *u L∞* |
| :---- | :---- | :---- |
| 0.0800000 0.0400000 0.0200000 0.0100000 0.0050000  | 0.0000000 0.9770784 1.1101687 1.2154725 1.5051878  | 0.00000001.0619210 1.0970910 1.0156016 1.3688708 |

T~~ABLE~~ 4.4  
*Table of EOC coefficients for the composite scheme.* 

*h* is presented in figures at time *t* \= 0\. The initial velocities *u*,*v* were set to zero. The solution was computed at the time interval \[0*,* 0*.*14\] using the space domain \[*−*2*,* 2\] *×* \[*−*2*,* 2\]. We use transmissive and reflective boundary conditions at points *Ub* on the boundary. The transmissive boundary condition are defined as *Ub* \= \[*h, uh, vh*\]*T*and the reflective ones are defined as *Ub* \= \[*h,* 0*,* 0\]*T*, where *h, u, v* are values at the inner neighbour of each boundary point.  
30 S. Brand Shallow water, t=0.0 

0.8 

0.6 

0.42 

Shallow water, t=0.02 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.04 

0.8 

0.6 

0.42   
Shallow water, t=0.06 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.08 

0.8 

0.6 

0.42   
Shallow water, t=0.10 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.12 

0.8 

0.6 

0.42   
Shallow water, t=0.14 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

FIG. 5.1. *Shallow water, transmissive boundary condition.*  
Parallel algorithm for numerical solution of problems in fluid dynamics 31 

Shallow water, t=0.0 

0.8 

0.6 

0.42   
Shallow water, t=0.02 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.04 

0.8 

0.6 

0.42   
Shallow water, t=0.06 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.08 

0.8 

0.6 

0.42   
Shallow water, t=0.10 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.12 

0.8 

0.6 

0.42   
Shallow water, t=0.14 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

FIG. 5.2. *Shallow water, transmissive boundary condition.*  
32 S. Brand Shallow water, t=0.0 

0.8 

0.6 

0.42 

Shallow water, t=0.02 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.04 

0.8 

0.6 

0.42   
Shallow water, t=0.06 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.08 

0.8 

0.6 

0.42   
Shallow water, t=0.10 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

Shallow water, t=0.12 

0.8 

0.6 

0.42   
Shallow water, t=0.14 

0.8 

0.6 

0.42 

\-1 \-2 

0 

1   
\-2 2 

0   
\-1 

1 

\-1 \-2 

0   
1 

\-2 2 

0   
\-1   
1 

FIG. 5.3. *Shallow water, reflective boundary condition.*  
Parallel algorithm for numerical solution of problems in fluid dynamics 33 

FIG. 6.1. *Used splitting of the space domain.* 

6\. Parallelization of numerical algorithms. The main purpose of our work was to compare the efficiency of parallel algorithms for numerical solution of the shallow water equation on systems with shared memory. Parallel computation is based on concepts for data exchange, shared and distributed memory. Shared memory means that all data are saved in the memory that can be accessed by all CPUs. This concept is used by OpenMP API (see \[2\]). Distributed memory means that in a multiprocessor system each processor has its own memory and after the processing of required computation tasks the data have to be reassembled. This concept is used by Message Passing Interface (MPI) (see \[1\]) and is applied similarly to \[8\]. 

All the presented numerical schemes are explicit. This enables to split the space domain into *p* parts and solve each part on one processing unit. The OpenMP shares the memory, therefore, after computing values at time *n*, all CPUs start computing the *n* \+ 1 time level. However, with distributed memory systems, it is necessary to send values from the boundary of each part to all neighbors after completing every time step computation. 

The domain splitting is illustrated in Figure 6.1, where the square represents the border of the space domain and dash lines represent borders of parts used for the computation on each processing unit. It can be seen that each middle part has only two neighbors. This simplifies communication needs for parallel computation. One can admit that this splitting is not minimizing the number of data sent. For example, if we use the domain splitting with the square parts, then we will need to send only two thirds of the data needed for this case. But this splitting is independent of the number of processing units, which allow us to measure the efficiency with the same program on any number of processors. 

7\. Results of the efficiency measurement. The efficiency measurement was performed for the shallow water equation (2.4) on the following parallel systems. The first three com puters are shared memory parallel systems and the last two computers are distributed parallel systems. 

*•* Hewlett Packard C8000, 4xCPU HP PA-RISC \- 1GHz, 12GB RAM *•* SGI ALTIX 3700, 24xCPU Intel Itanium II \- 1,3 GHz, 64GB RAM *•* Linux PC, 2xCPU Intel Pentium II \- 700 MHz, 1GB RAM 

*•* Linux grid, 8xCPU Intel Pentium 4 \- 2.4 GHz, 512MB RAM 

*•* Linux grid, 12xCPU Intel Pentium 4 \- 2.4 GHz, 512MB RAM 

The results of efficiency measurement are presented in Tables 7.1 and 7.1. The tables have the following structure. The first column contains the grid dimension. The second  
34 S. Brand 

column contains the time of sequence program in seconds, this means the time of the compu tation made by only one processing unit. The remaining columns contain the time of parallel program and the efficiency of this program in the brackets. In the header of these columns there is specified how many MPI and OpenMP processes were used in the computation. The time duration of each computation was measured by the C *gettimeof day*() function as a difference between the start and the end time. The times listed in here are the times needed for the computation only. This means the times needed for the value initialization and result saving is excluded. The efficiency is calculated from the following formula: 

efficiency \= sequence time   
parallel time *×* number of processors (7.1) 

| grid  |  | OMP=4,MPI=1  | OMP=2,MPI=2  | OMP=1,MPI=4 |
| ----- | ----- | ----- | ----- | ----- |
| 100*×*100  200*×*200  300*×*300  400*×*400  500*×*500  600*×*600  700*×*700  800*×*800  900*×*900  1000*×*1000 1100*×*1100 1200*×*1200 1300*×*1300 1400*×*1400 1500*×*1500 1600*×*1600 1700*×*1700 1800*×*1800 1900*×*1900 2000*×*2000  | 0.660  5.705  19.436  46.568  89.915  154.264  248.072  369.939  533.809  726.617  960.135  1250.901  1600.657  1975.505  2432.951  2918.703  3562.641  4169.705  4900.167  5668.034  | 0.661(0.250)  1.704(0.837)  5.445(0.892)  12.432(0.936)  25.429(0.884)  43.074(0.895)  66.403(0.934)  97.479(0.949)  138.156(0.966) 194.368(0.935) 252.363(0.951) 327.944(0.954) 432.835(0.925) 530.362(0.931) 639.059(0.952) 762.497(0.957) 943.562(0.944)  1142.250(0.913) 1289.412(0.950) 1513.491(0.936)  | 0.249(0.664)  1.833(0.778)  5.525(0.879)  12.564(0.927)  24.972(0.900)  42.324(0.911)  66.086(0.938)  95.952(0.964)  139.469(0.957) 190.263(0.955) 248.010(0.968) 322.275(0.970) 429.117(0.933) 523.026(0.944) 635.187(0.958) 754.847(0.967) 931.938(0.956)  1139.587(0.915) 1278.466(0.958) 1495.630(0.947)  | 0.262(0.631) 1.540(0.926)  5.283(0.920)  12.411(0.938)  23.695(0.949)  40.290(0.957)  63.709(0.973)  94.033(0.984)  134.933(0.989) 183.260(0.991) 244.049(0.984) 319.292(0.979) 408.119(0.981) 506.207(0.976) 626.386(0.971) 753.314(0.969) 914.391(0.974)  1078.108(0.967) 1269.950(0.965) 1462.494(0.969) |

T~~ABLE~~ 7.1  
*Time and efficiency of parallel program using composite scheme computed on sugar.fjfi.cvut.cz.* 

6000 5000 4000 3000 2000 1000 0   
computation time 

OMP=1,MPI=1   
OMP=4,MPI=1   
OMP=2,MPI=2   
OMP=1,MPI=4 

2000150010005000 

2 

1.5 1 

0.5 0  
efficiency 

OMP=1,MPI=1   
OMP=4,MPI=1   
OMP=2,MPI=2   
OMP=1,MPI=4 

2000150010005000   
Parallel algorithm for numerical solution of problems in fluid dynamics 35 

| grid  |  | OMP=2,MPI=1  | OMP=1,MPI=2 |
| ----- | ----- | ----- | ----- |
| 100*×*100  200*×*200  300*×*300  400*×*400  500*×*500  600*×*600  700*×*700  800*×*800  900*×*900  1000*×*1000 1100*×*1100 1200*×*1200 1300*×*1300 1400*×*1400 1500*×*1500 1600*×*1600 1700*×*1700 1800*×*1800 1900*×*1900 2000*×*2000  | 0.660  5.705  19.436  46.568  89.915  154.264  248.072  369.939  533.809  726.617  960.135  1250.901  1600.657  1975.505  2432.951  2918.703  3562.641  4169.705  4900.167  5668.034  | 0.407(0.812)  2.979(0.957)  10.542(0.922)  23.888(0.975)  47.564(0.945)  78.828(0.978)  128.731(0.964) 185.831(0.995) 274.079(0.974) 364.812(0.996) 490.989(0.978) 629.919(0.993) 815.476(0.981) 994.999(0.993)  1240.065(0.981) 1470.887(0.992) 1809.891(0.984) 2101.875(0.992) 2491.581(0.983) 2858.638(0.991)  | 0.350(0.943) 3.018(0.945)  10.168(0.956)  24.081(0.967)  46.331(0.970)  78.992(0.976)  125.554(0.988) 187.054(0.989) 267.076(0.999) 365.749(0.993) 484.739(0.990) 632.366(0.989) 809.857(0.988)  1002.577(0.985) 1238.091(0.983) 1484.173(0.983) 1799.128(0.990) 2113.957(0.986) 2485.624(0.986) 2891.638(0.980) |

T~~ABLE~~ 7.2  
*Time and efficiency of parallel program using composite scheme computed on sugar.fjfi.cvut.cz.* 

6000 5000 4000 3000 2000 1000 0   
computation time 

OMP=1,MPI=1   
OMP=2,MPI=1   
OMP=1,MPI=2 

2000150010005000 

2 

1.5 1 

0.5 0   
efficiency 

OMP=1,MPI=1   
OMP=2,MPI=1   
OMP=1,MPI=2 

2000150010005000 

8\. Conclusion. The efficiency results show that both OpenMP and MPI are suitable for paralleling programs computing numerical solution of the shallow water equation because its values for the grid that was large enough were never less than 85% during the testing. For a small grid there the communication time is higher than the processing time and therefore the efficiency is too small. One important result is that efficiency of the program using MPI was higher than efficiency of the program using OpenMP. This is due to the ideal conditions over the measuring. Any load of one processing unit will enlarge the computational time of the program using MPI while the program using OpenMP will be minimally affected. OpenMP contains built-in tools for load balancing that can be easily used even if the programmer did not originally intend to use it. The MPI programmer has to write his own load balancing which will enlarge the programming time.  
36 S. Brand 

Acknowledgment. The work was partly supported by the project No. LC06052 ”Jindˇrich Necas Center for Mathematical Modelling” of the Ministry of Education, Youth and Sports ˇ of the Czech Republic. 

REFERENCES 

\[1\] *MPI: A Message-Passing Interface Standard*. http://www.mpi-forum.org/docs/mpi-11-html/mpi-report.html. \[2\] *OpenMP Application Program Interface*. http://www.openmp.org/drupal/mp-documents/spec25.pdf. \[3\] D. L. GEORGE, *Numerical Approximation of the Nonlinear Shallow Water Eguations with Topography and* 

*Dry Beds: A Godunov-Type Scheme*. University of Washington, Washington (2004). Master Theses. \[4\] R. J. LEVEQUE, *Numerical Methods for Conservation Laws*. Birkhauser Verlag, Basel (1990). \[5\] R. J. LEVEQUE, *Nonlinear Conservation Laws and Finite Volume Methods for Astrophysical Fluid Flow*. Springer-Verlag, University of Washington (1998). 

\[6\] R. LISKA, AND B. WENDROFF, *Composite Schemes for Conservation Laws*. SIAM J. Numer. Anal. 35, no. 6 (1998), 2250–2271, . 

\[7\] R. LISKA, AND B. WENDROFF, *2D Shallow Water Equations by Composite Schemes*. Int. J. Numerical Methods in Fluids 30 (1999), 461–479.   
\[8\] M. Sˇ ENKY´ Rˇ, J. MIKYSKA ˇ , AND M. BENESˇ, Application of Parallel Computing Techniques for Problems of Degenerate Diffusion. In *Numerical Mathematics and Advanced Applications, ENUMATH 2003 (peer reviewed proceedings)*, pages 756–766, Springer Verlag, Berlin, 2004\. eds. M. Feistauer, VDolejsˇ´ı, P. Knobloch, K. Najzar, ISBN 3-540-21460-7.