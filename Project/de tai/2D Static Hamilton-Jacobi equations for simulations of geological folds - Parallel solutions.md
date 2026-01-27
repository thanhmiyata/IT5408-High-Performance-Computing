Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10   
http://www.mathematicsinindustry.com/content/4/1/10 

**R E S E A R C H Open Access** 

Parallel solutions of static Hamilton-Jacobi equations for simulations of geological folds 

Tor Gillberg1,2\*, Are Magnus Bruaset1, Øyvind Hjelle2 and Mohammed Sourouri1 

\*Correspondence: torgi@simula.no 1Simula Research Laboratory, Martin Linges vei, Oslo, Norway   
2Kalkulo AS, Martin Linges vei, Oslo, Norway   
**Abstract** 

Two new algorithms for numerical solution of static Hamilton-Jacobi equations are presented. These algorithms are designed to work efficiently on different parallel computing architectures, and numerical results for multicore CPU and GPU implementations are reported and discussed. The numerical experiments show that the proposed solution strategies scale well with the computational power of the hardware. The performance of the new methods are investigate for tow types of static Hamilton-Jacobi formulations are investigated, the isotropic eikonal equation and an anisotropic formulation used to simulate different types of geological folding. Simulations of geological folding is a key component in an industrial software used in oil and gas exploration. In particular, our experiments indicate that the new algorithms would be capable of accelerating an existing industrial simulator substantially. Direct comparison with the current industry code shows that computing times can be reduced from several minutes to a few seconds. The new methods are now being migrated to the industrial software. 

**Keywords:** parallel computing; general purpose GPU computing (GPGPU); multicore CPU; eikonal equation; static Hamilton-Jacobi equations; front propagation; geological folding; generalised distance computing; CUDA; OpenMP 

**1 Static Hamilton-Jacobi equations** 

This paper concerns methods for simulating monotonically propagating fronts modelled by static Hamilton-Jacobi equations. A general formulation of such equations is 

*H* **x**,∇*T*(**x**) \= , 

*T*(**x**) \= *g*(**x**) ∀**x** ∈  , 

() 

where *g*(**x**) represents initially given values at points in the set  , and the Hamiltonian *H* is here assumed Lipschitz continuous and convex in ∇*T*. The unknown entity, *T*, can be thought of as a generalised distance field, as well as the time of arrival of a monotonically expanding front. Usually, *g*(**x**) \=  which means that   represents the object from which the generalised distance should be computed. Static Hamilton-Jacobi equations are known to often have multiple-valued solutions, in that a point in the domain may have several values simultaneously. The viscosity solution of the static Hamilton-Jacobi equations is also the minimal distance value, that is the first time of arrival \[\]. In some applications, the multiple-valued solutions can be important \[\]. However, in many applications it is 

©2014 Gillberg et al.; licensee Springer. This is an Open Access article distributed under the terms of the Creative Commons Attribution License (http://creativecommons.org/licenses/by/2.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original work is properly cited.  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 2 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

enough to solve for the viscosity solution of the equation. In this paper we solve solely for 

the viscosity solution. 

**1.1 A mathematical framework for propagating fronts** 

The solution *T* of () can be thought of as the time of arrival of a wavefront that propagates 

through a domain. More generally, the *t*\-isosurface of *T* shows the position of the front 

at time *t*. The *t*\-isosurface can also be regarded as the generalised *t*\-distance from  , that 

is, solution values are increasing away from  . In a front propagation setting, equation () 

is often formulated as *V* ∇*T*  \= , where the scalar *V* is the speed at which the front is 

moving in the direction normal to the front. At every point, *V* is greater or equal to , 

with equality implying that the front never reaches the point. If the speed is independent 

of the direction the front is moving in, the front propagation is said to be isotropic. The 

isotropic eikonal equation is important in many applications and is formulated as 

*F*(**x**)  ∇*T*(**x**)   \= . () 

For the special case of *F* ≡  and *g*(*x*) ≡ , the viscosity solution of the eikonal equation is 

the minimal Euclidean distance to  . 

The problem is anisotropic if the rate of increase of *T* values changes with the direction 

it is measured in. There is a clear causality interpretation of a propagating front since the 

current motion does not depend on future events. In other words, smaller values of *T* 

are independent of larger *T* values, which imply that the discrete solution representing *T* 

should be computed in increasing order. In the isotropic eikonal case this causality obser 

vation can be directly translated to discrete nodal values. However, for the more general 

anisotropic case a similar discrete interpretation is not valid \[\], and anisotropic problems 

are therefore more complicated to solve. The algorithms developed in this paper have been 

applied to both isotropic eikonal formulations and to anisotropic problems of the form (). 

Most formulations of static Hamilton-Jacobi equations are nonlinear, and the compu 

tations needed for updating the solution in a grid node are expensive in terms of both 

floating-point operations and logical branching. The amount of computations needed can 

be decreased by updating nodes in different orderings. 

In the following section we give a brief overview of solution algorithms. For a more 

complete algorithmic overview the reader is to \[\] and references therein. 

**1.2 Applications** 

Many physical phenomena can be described by a propagating front. A front can for ex 

ample describe an interface between different objects or fluids (multi-phase flow) \[\], a 

shock wave \[\], or the arrival of a wave \[\]. In this work, the front is assumed mono 

tonically expanding, like a wildfire spreading only to unburned grounds. Monotonically 

expanding fronts can be entirely described by the time of arrival of the front to all points 

in the domain. Applications requiring fast simulations of monotonic front propagations 

include medical tomography \[\], simulation of cardiac activation times \[\], segmentation 

of images \[\], and seismic wave propagation \[\]. Shorter computing times may result in 

faster and more accurate models since one can afford to better explore parameter ranges, 

and thereby the model uncertainty. As a result, faster methods can provide further insight 

into the underlying physical problems. In this paper, we focus on the development of new  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 3 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

efficient algorithms. For a proper introduction to front propagation and a range of appli 

cations, see \[\]. The application initiating the research presented in this paper is from a 

software project developed by Statoil ASA and Kalkulo, and is presented in the following 

paragraph. 

*The compound earth simulator* 

Structural restoration and reconstruction is a methodology for building and improving 

structural models of subsurface geology. The Compound Earth Simulator (CES) is a soft 

ware in which large three dimensional models can be transformed between different struc 

tural states \[\]. Generalised distance fields are needed in the transformations of the geol 

ogy and to deposit the layers. The generalised distance fields are given from the following 

static Hamilton-Jacobi equation, 

*F*  ∇*T*(**x**)   \+ *ψ* **a** · ∇*T*(**x**) \= , **x** ∈  , *T*(**x**) \= *t* ∀**x** ∈  . 

() 

Here,   is the computational domain,   is the surface of a geological layer (horizon), and **a** is a unit vector. Geological volumes are transformed repeatedly in the workflow \[\]. For the software to be user friendly, the distance fields must be computed rapidly while still being accurate. 

In this paper, software is said to be interactive if the delay due to computations is of the order of a few seconds at the most. Some of the available algorithms are fast enough for interactive use in two spatial dimensions. This is not the case in three dimensions, since the number of nodes and the nodal update cost increase significantly. A case study of geological folding based on structures extracted from offshore seismic data, is included in the numerical experiments reported in this paper. The methods presented in this study enables interactive high-resolution restorations in three dimensions. 

*The ray equations* 

Characteristic curves are paths along which ‘particles’ on the front are pushed forward, as well as lines along which distances are measured. Let **x**(*s*)=(*x*, *y*, *z*)(*s*) denote a char acteristic curve parameterised by the scalar *s*. Characteristic curves are often referred to simply as characteristics. The solution *T* is strictly increasing along characteristics that are perpendicular to the moving front. The characteristic equations can be formulated on basis of () as 

*d***x**   
*ds* \= ∇**p***H*(**p**, **x**), () *d***p**   
*ds* \= ∇**x***H*(**p**, **x**), () 

where **p** \= ∇*T*, ∇**x** \= (*∂*/*∂x*, *∂*/*∂y*, *∂*/*∂z*)*T* , and ∇**p** \= (*∂*/*∂Tx*, *∂*/*∂Ty*, *∂*/*∂Tz*)*T* \[, \]. For the Hamiltonians () and (), the respective characteristic curves **x***e*(*s*) and **x***f* (*s*) are given by 

**x***e*(*s*) \= **x***e*() \+ *sF***n**, () **x***f* (*s*) \= **x***f* () \+ *s*(*F***n** \+ *ψ***a**), ()  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 4 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

where **n** \= ∇*T*/ ∇*T*  is the unit normal vector to the propagating front, and *x*() marks 

the starting location of the curve. The velocity parameters *F* and *ψ***a** are assumed to be 

constant in a small neighbourhood, which results in locally linear characteristic curves 

along which **p** is constant. On a discrete grid, the parameters are commonly assumed to 

be constant within the stencil area. An investigation of the solution value along a charac 

teristic curve shows that 

*dT*   
*ds* \= **p** ·*d***x***ds* \= **p** · ∇**p***H* \=  for () and () ()   
 ⇒ *T* **x**(*s*) \= *T* **x**() \+ *s*. () 

That is, *T* is increasing linearly in *s* along characteristics. Because of this one-dimensional 

dependency, solution values can be extrapolated with high accuracy along characteris 

tic curves. An understanding of characteristics is therefore important when creating nu 

merical stencils, as well as designing new algorithms for front propagation problems. The 

solutions *Te* and *Tf* for equations () and () change explicitly according to 

*Te*(*s***n**) \= *Te*() \+ *sF* , () 


*Tf*   
*sF***n** \+ **a**  *F***n** \+ **a**    
    
\= *Tf* () \+ *s* 

 *F***n** \+ **a**  () 

along characteristics. In particular, *Te* increases at a constant rate since the front moves with a speed that is independent of the normal and its direction. This is not the case for *Tf* , which increases with a rate that depends on the direction in which the front moves. From these formulations the increased complexity for the anisotropic formulation is obvious. The anisotropic dependencies in the fold equation complicate the computations, and the computational cost is significantly higher than for the isotropic problem. 

**1.3 Numerical solution methods** 

The solution can be built by computing the local solution along a set of characteristic curves, and thereafter interpolating the travel time between the curves. This approach is traditionally referred to as the method of characteristics, or ray tracing. In seismology, the term ‘rays’ is often used in lieu of characteristics. The methods developed in this paper assumes a fixed grid of nodal values. Grid-based methods use stencils that are strongly related to ray-tracing methods. Instead of tracing values along specific characteristics, the grid-based tracking methods update nodal values by interpolating solution values be tween neighbouring nodes organised in stencils. Commonly, solution values are interpo lated linearly along edge-connected nodes in the grid \[\], implying a planar modelling of the front. In Appendix , we present an efficient stencil formulation for the eikonal equa tion. A smaller value corresponds to an earlier arrival time, and thus a better estimate of the sought viscosity solution. The solution at every point is therefore strictly decreasing. A too small value will not be enlarged, but instead it generates an error that will spread throughout the domain. It is therefore important that new estimates are based on upwind information, that is, the characteristic originates from values that are upwind of the node being updated and therefore already passed by the front.  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 5 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

![][image1]![][image2]  
(a) Tracking order. (b) Parallel marching sweep. 

**Figure 1 Illustration of update order of algorithms, where circles are nodes on the grid.** (**a**) Tracking   
order. Red nodes belong to the narrow band, describing the current front position. One node at a time is   
moved from the narrow band to the collection of already passed nodes (in blue). Nodes that have yet to be   
assigned a first value are empty circles. (**b**) The four sub-sweeps of the parallel marching method. Red nodes   
are being computed using values of nodes connected by arrows, which points in the direction of the   
sub-sweep. Nodes are coloured darker after each time their value has been updated. 

*Tracking methods* 

Methods that compute nodes in the same order as they are reached by the front, are re 

ferred to as front *tracking* methods. Tracking algorithms classify a node as passed by the 

front only once, and are often referred to as one-pass algorithms. Therefore, these al 

gorithms compute solution estimates only a few times per unknown node. The current 

position of the front is described by a set of nodes known as the *narrow band*. Figure (a) 

is an illustration of the classification of nodes in tracking methods. To know which value 

is the to be considered passed by front, the narrow band nodes are stored in a min-heap 

data structure that is kept sorted \[, \]. With *N* being the total number of nodes, the 

worst case computational cost of the algorithm is *O*(*N* log (*N*)), where the log (*N*) factor is 

due to the min-heap sorting. It is important to note that a direct application of the causal 

ity principle to discrete nodal values does not hold for anisotropic problems \[, \], and 

the original tracking method sometimes fails to converge \[\]. An algorithmic extension 

for tracking solvers, is to use adaptive (also called dynamic) stencil shapes \[, \]. When 

updating a node, the stencil shape then varies depending on the dynamics of the front, 

and nodes covering a wider area must contribute to the update step. The footprint of the 

stencil is often determined on basis of the anisotropy coefficient \[\], which is a measure 

of the degree of anisotropy in the chosen problem. 

The need for a sorted data structure makes tracking methods conceptually sequential. 

However, there are some parallel implementations of tracking-like methods. The group 

marching method \[\] classifies several nodes as passed by the front in every iteration. 

Since several nodes are passed simultaneously, the method allows a certain degree of par 

allel processing. However, a parallel implementation is not straightforward since all neigh 

bours should be updated according to two different orderings to ensure convergence. The 

entire domain can also be decomposed into subdomains in which separate processors up 

date nodes with tracking methods \[, \]. Subdomains must then be padded with extra  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 6 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

layers of nodes, referred to as *ghost nodes* or *ghost points*. In order for rays to travel be 

tween subdomains, the resulting algorithms are not strictly of tracking type, since already 

computed subdomains can be recomputed several times. 

*Iterative methods* 

Sweeping methods are iterative methods based on the observation that the minimal dis 

tance follows along a unique direction. By sweeping through the grid in one direction at 

a time, the minimal distances in that direction is computed for all nodes \[\]. Computed 

nodes are used shortly after they have been assigned new values. The characteristic curves 

are extended in the same direction as the sweep moves through the domain. Moreover, 

sweeping methods are almost directly applicable to anisotropic problems \[\], and can 

be adapted to parallel architectures with good performance on multicore systems with a 

domain-division approach \[, \]. Still, the parallelism of traditional sweeping methods 

is coarse grained. This is caused by neighbouring nodes in a sweep depending on each 

other’s values. Nodes on ‘slices’ of the domain do not share dependencies, and can there 

fore be computed in parallel \[\]. Slices are small in corners of the domain, and parallel 

processing is then not beneficial. Another disadvantage is that memory accesses of nodes 

on a slice are incoherent, which impedes the performance of GPU implementations \[\]. 

The parallel marching methods also sweep through the domain, but in slightly modi 

fied directions and with alternative stencil shapes \[, \]. The alternative stencil shapes 

remove dependencies between neighbouring nodes, so that they can be updated simul 

taneously. Figure (b) give an illustration of the sweeping process of the two dimensional 

parallel marching method. A full sweep of the grid consists of four sub-sweeps, during 

which the front is propagated up, right, down and in the left direction. In two spatial 

dimensions, lines of nodes can be updated simultaneously, and in three dimensions an 

entire layer of nodes can be updated in parallel \[, \]. The method scales well when im 

plemented on GPUs. Sweeping methods are faster than tracking methods when the do 

main structure and the velocity formulation are simple. For complicated problems, many 

sweeps are needed for convergence and the methods become slow \[, \]. Similarly, when 

anisotropy is dominant, more iterations are often needed for convergence \[\]. Traditional 

sweeping methods need to sweep the domain in *k* directions in *k* dimensions \[\] to span 

the solution space. Interestingly, only *k* sweeps are needed for extensions of the parallel 

marching method to higher dimensions. However, in practice more full sweeps of the do 

main might be needed by the parallel marching method than by the traditional sweeping 

methods \[\]. 

*Label correcting methods* 

There are algorithms that combine iterative and tracking methods by tracing a front with 

a less strict ordering. A common name for such methods is *label correcting*, emanating 

from shortest path methods in graph networks \[, \]. Tracking methods are referred to 

as *label setting*, since nodes are labelled as ‘passed’ only once. Since label correcting meth 

ods use a less strict ordering than tracking methods, they can sometimes be implemented 

on parallel architectures. For instance, the massive marching method \[\] updates an en 

tire list of nodes in parallel on a multicore CPU. If a node receives a new value during the 

update, all direct neighbours with larger solution values are collected in a new list. The 

procedure is repeated until there are no nodes to be added to a list. Massive marching  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 7 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

is very similar to the fast iterative method \[\] (FIM), which has also been used to solve 

problems with anisotropy \[\]. GPU FIM have been used in several applications \[, \], 

and can be extended to compute geodesic distances on triangulated surfaces \[\]. Large 

velocity contrasts may cause the list to grow long and nodes to be recomputed many times, 

thereby increasing the computing time significantly \[\]. The amount of needed compu 

tations is reduced if nodes are updated in an order more similar to the actual position of 

the expanding front. The two-queue method \[\] and semi-ordered FIM (SOFI) \[\] are 

more uniform in their performance since they enforce a stricter ordering of the updates. 

In these methods, only selected neighbours are directly added to the new list, leaving for 

later computations the ones too far ahead of the front. As a result, the lists in SOFI and 

the two-queue method follow the isocurves of the solution more closely than FIM. 

**2 Developing fast parallel solvers** 

For a front propagation method to perform efficiently, the algorithm must make use of 

the causality information embedded in the problem. Brute force iterative and sweeping 

approaches are too costly in terms of computations, and methods that are strictly track 

ing the front are too sequential. For instance, the parallel marching algorithms work well 

on GPUs, but their performance on multi-core devices is not impressive when compared 

to the algorithms presented in this paper. Since the entire grid is updated repeatedly, such 

sweeping-like methods cannot perform well enough due to the high amount of compu 

tations. Tracking methods are stable, and each node is recomputed only a few times \[\]. 

Since tracking methods are conceptually sequential, the computing time increases no 

ticeably when the grid size increases. Extensions for problems with anisotropy are com 

plicated and sometimes not applicable \[\]. In contrast, the more general class of label cor 

recting methods can solve anisotropic problems \[\] and can sometimes be implemented 

on parallel computing architectures \[, \]. Such methods are therefore particularly in 

teresting to investigate further. 

**2.1 Algorithmic inspiration** 

Several recent contributions to the field have inspired the design of the algorithms pre 

sented in this paper. In this section we discuss the relation between our algorithms and 

other existing methods. 

*Two-scale methods* 

For problems where the characteristics are straight lines, sweeping methods converge 

quickly and are often faster than tracking methods. However, problems with highly non 

linear characteristic curves often require many sweeps to reach convergence, and ordered 

methods then perform better. For most problem formulations, the characteristics bend 

smoothly and can be approximated well as piecewise linear.With these properties in mind, 

Chacon et al. \[\] design several two-scale algorithms for isotropic problems. In the two 

scale methods, the domain is split into subdomains and an ordering of the subdomains 

is created, for instance with a FMM solver acting at the subdomain level. Every subdo 

main is thereafter computed with a sweeping method, one subdomain at a time in the 

determined order. Alternatively, subdomains can enter a heap repeatedly in the *heap-cell* 

*method* (HCM). The local sweeping method converges fast since characteristic curves do 

not bend much within a neighbourhood, such as a subdomain. Recently, the same authors  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 8 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

introduced a parallel two-scale method with the *parallel Heap-Cell Method* (pHCM) for 

eikonal equations \[\]. The pHCM targets specifically multicore architectures, and show 

impressive speedups when compared with other algorithms. Every thread has its own heap 

with non-overlapping subdomains awaiting computations. In this paper we consider more 

general anisotropic front propagations. 

*Domain decomposition* 

Two subdomains can be updated independently of each other as long as they do not share 

a boundary. On shared boundaries there is a risk of read-write memory interferences \[\]. 

To avoid any collision between read and write operations, all subdomains are padded with 

a layer of ghost nodes \[\]. Ghost nodes are local copies of nodes belonging to subdomains 

that are immediate neighbours. After a subdomain has been computed, ghost nodes in 

surrounding subdomains must be updated with the most recently computed values. This is 

done in a *synchronisation* step that follows the computation of a subdomain. The boundary 

condition for a subdomain is then based on the values that are synchronised from the 

neighbouring subdomains. Subdomains surrounded by ghost nodes are further discussed 

in Appendix . Any algorithm can be used to compute new values within a subdomain. 

*Locks* 

The concept of locking is a methodology for reducing the number of unnecessary compu 

tations, which was introduced at the nodal level in the locked-sweeping method \[\]. The 

solution value of a node needs to be updated only if values on neighbouring nodes have 

changed since the last update. Therefore, a node can be locked for computation right after 

an update, and be kept locked until relevant changes are observed in its vicinity. A lock 

can be considered a label that can be reset several times. Consequently, algorithms using 

locks are of the label correcting type. We have extended the idea of locks to the subdomain 

level. In our algorithms, we use a compute lock (`CL`) to prevent unnecessary computations 

of subdomains. After computing all values in a subdomain, it is locked for computation. 

This compute lock is unlocked whenever the subdomain receives a new value in the syn 

chronisation step. Adopting the terminology used in FIM, we refer to a subdomain with 

an open compute lock as *active*. 

**2.2 Two new algorithms** 

In this paper, we propose two new algorithms that process sets of subdomains in parallel. 

We refer to the sets of subdomains as *schedules*. The order in which subdomains are com 

puted depends on which subdomains are put into schedules. A schedule is represented as a 

list of indices to subdomains. The new algorithms differ mainly in the way they build their 

respective schedules. All subdomains in the schedule are computed and synchronised re 

peatedly until most subdomains have been locked for computations, and a new schedule 

is built. In the development of these algorithms also a third algorithm was investigated, 

called the *two-scale parallel marching method* \[\]. Since this method is outperformed by 

the list-based methods, a presentation is omitted from this paper. 

In the reported implementations, we compute new nodal values in a subdomain with the 

*three-dimensional parallel marching method* (D PMM) \[\], of which details are given 

in Appendix B.. Since D PMM is a parallel method, the resulting algorithm has two lev 

els of parallelism. Hardware support for multiple layers of parallelism exists for instance  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 9 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

in GPUs. More advanced subdomain solvers (FMM, FIM or SOFI) may increase perfor 

mance, especially when implemented on devices with lesser processing capabilities. How 

ever, such investigations are beyond the scope of this paper. The remainder of this section 

presents the two new algorithms in detail. 

*Algorithmic framework* 

The computed values *Ti*,*j*,*k* approximating *T* in the nodal positions (*i*, *j*, *k*) are stored in 

an array corresponding to a rectangular grid in three spatial dimensions. The domain is 

decomposed in (*sx*,*sy*,*sz*) subdomains in the *x*, *y*, and *z*, directions, respectively. Every sub 

domain consists of (*bx*, *by*, *bz*) nodes, making the number of nodes in the *x* direction of 

the rectangular domain *bxsx*. All subdomains are surrounded by a layer of ghost nodes on 

either sides. The ghost nodes are stored physically in memory, resulting in a total stor 

age of (*bx* \+ )*sx* nodal values along the *x*\-axis. Similarly there are (*by* \+ )*sy* and (*bz* \+ )*sz* 

nodes along the two other axes when ghost nodes are included. The grid may have a non 

uniform spacing defined by the parameters *dx*, *dy*, and *dz*, but the same spacing values 

apply to all subdomains. When algorithmic details are discussed on a nodal basis, nodes 

are indexed with (*i*, *j*, *k*) while subdomains are indexed with (*ii*, *jj*, *kk*). For increased read 

ability, subdomains are sometimes indexed by the triple `sd` \= (*ii*, *jj*, *kk*). The compute lock 

`CL` is implemented as an array with one binary state value for each subdomain, `Open` or 

`Locked`. 

Similar to most front propagation methods, the proposed algorithms assume mono 

tone convergence of the solution. All nodes are first assumed never to be reached by the 

front, *Ti*,*j*,*k* \= ∞. A set of nodes with initial values is assumed known initially. These values 

constitute the discretised version of the boundary condition in (). The subdomains that 

receive at least one initial value are unlocked for computations during the initialisation 

phase. As new *T* values are computed, they are only accepted if they are smaller than the 

previous value, thus assuring monotone convergence from above. The solution is bounded 

from below by the smallest initialised value since the solution *T* is increasing away from 

 . Note that solution values will not be increased during the iterations once accepted, 

and too small values should therefore never be accepted. Too small approximations can 

be avoided in upwind stencils if the used solution values lie upwind of the node being up 

dated. Further details on the construction of stencils for front propagation problems are 

discussed in Appendix , where an efficient eikonal stencil is presented. 

*Generic solver* 

The proposed algorithms have a top level structure given by the template in Algorithm . 

As discussed below, each algorithm have individual specifications of the conditions gov 

erning the loops marked *A* and *B* . After initiating the problem, a schedule of computing 

tasks is created. In this context, a task corresponds to the computations needed for a single 

subdomain. The scheduled tasks are computed and the ghost nodes updated repeatedly 

through synchronisations with neighbouring subdomains, before a new schedule is cre 

ated. 

The compute and synchronisation steps are similar for both algorithms while the actual 

construction of the schedule differs. After a subdomain `sd` has been computed, the com 

pute lock is closed, `CL`(`sd`) ← `Locked`. In ComputeSchedule(), all active scheduled 

subdomains are computed in parallel as shown in Algorithm .  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 10 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Algorithm**  Solver() 

**Algorithm**  ComputeSchedule(`Schedule`) 

**Algorithm**  SyncFromSchedule(`Schedule`) 

The ComputeSubdomain() function updates the solution for all nodes and ghost 

nodes in the subdomain `sd`. Because of the presence of the ghost nodes, the compu 

tations rely only on nodal values stored locally for the subdomain, thus eliminating the 

need for subdomain communication during the computations. After computing new val 

ues for the active subdomains in the schedule, nodal values must be communicated be 

tween adjacent subdomains to assure that all ghost nodes have correct values. This is done 

in SyncFromSchedule(), shown in Algorithm . 

To avoid any memory interference, the nodal values of neighbouring subdomains are up 

dated in parallel in one direction at a time. First, the ghost nodes belonging to subdomains 

with larger *i* values are updated in parallel, and thereafter ghost nodes in subdomains with 

larger*j*, *k*, and smaller*i*, *j*, *k* indices are updated. Note that no computations are performed 

during the synchronisation, rather the nodal values are set to match their copies in the 

other subdomain. A nodal value is only changed if the value of its counter-node is smaller.  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 11 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Algorithm**  BuildScheduleLAS(`L`,`CL`) 

If a node in subdomain `sd` receives a new value, the subdomain is activated by opening 

the compute lock, `CL`(`sd`) ← `Open`. Further details of the Sync Values() functions are 

given in Appendix B.. The algorithms finish when all subdomains are locked for compu 

tations and the schedule is empty. 

*.. The list of active subdomains method* 

The list of active subdomains (LAS) algorithm schedules all active subdomains for com 

putation by storing their indices in a list. As the algorithm iterates through the grid, the 

list follows a front of solution estimates as it travels through the domain, similarly to the 

lists in the massive marching and FIM methods \[, \]. Subdomains that have already 

been computed may at a later stage be added to a new version of the list. Pseudocode for 

building the schedule (list) for the LAS method is shown in Algorithm . The returned 

variable `noSched` is the number of subdomains in the list `L`. 

After computing and synchronising the list of subdomains, the number of active sub 

domains left in the list are computed, here referred to as `noActive`. The schedule is 

computed again if the density of active subdomains is higher than a given fraction. Let 

`noSched` denote the length of the list. After some experimentation with *B* in Algo   
rithm , we have chosen to repeat the computations while  \< `noActive`   
~~`noSched`~~ . This thresh   
old gives good performance for most examples. Moreover, let `noC` denote the number 

of subdomains that can be computed simultaneously by the current computational plat 

form. If there are fewer active subdomains than what can be computed simultaneously, 

i.e. `noActive` \< `noC`, a new schedule is created. This optional condition to improve load 

balancing is included in the LAS pseudocode given in Algorithm , which replaces the 

generic Solver() in Algorithm . 

The LAS method has some similarities with FIM, but there are differences to notice: 

For FIM, only the GPU implementation computes subdomains, and no compute locks are 

used. In FIM, the scheduled tasks are computed only once before a new schedule is built, 

whereas in LAS the schedule is reused until at most (around %) scheduled subdomains 

are locked for computation. Any method can be used to compute a subdomain in LAS, 

whereas in FIM all nodes of a subdomain are computed iteratively until convergence. This 

brute force approach of FIM would likely reduce the performance on multicore systems. 

*.. The semi-ordered list of active subdomains method* 

In the LAS method, previously active subdomains may be reactivated at later stages of 

the algorithm. In some front propagation problems, this reactivation can force the list 

to grow long. A similar observation holds for FIM and SOFI, although SOFI performs  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 12 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Algorithm**  SolverLAS() 

more uniformly thanks to the enforced heuristic ordering of updates \[\], which makes 

the behaviour of the list more similar to the narrow band of tracking methods. The semi 

ordered list of active subdomains (SOLAS) method similarly schedules subdomains in lists 

that in a sense stay close in shape to isosurfaces of the solution. If a subdomain is too far 

ahead of the others, it is not added to the schedule. 

In order to measure the approximate location of a subdomain, each subdomain has an 

associated scalar value stored in a global array structure, here referred to as `SD`. The `SD` 

value for a subdomain is the smallest new *T* value that activates the subdomain, and is 

set during the synchronisation procedure. This minimal *T* value is similar to the ‘rollback 

state’ of the domain decomposition FMM \[\]. Details on handling the `SD` values in the 

synchronisation is given in Appendix B.. 

The `SD` value is an approximate position of the current front as it reaches a subdomain. 

Let `Av` denote the current average `SD` value of all active subdomains, while `oldAv` denotes 

the same average of the previous iteration. All subdomains with a `SD` value smaller than a 

cutoff value `cutT` are scheduled for computation. In our implementation the `cutT` value 

is set to be `Av` \+ .max (,`Av` – `oldAv`). Several methods of creating this cut-off value 

was tested, but the chosen . relaxation from the average value gives good performance 

for many problems. With the enforced ordering of subdomains, the list will mimic the 

behavior of the solution’s isosurfaces. 

If the list is made too short, some processing units may be idling since there are not 

enough computations to perform. Since the list is computed repeatedly, neighbouring sub 

domains may be reactivated repeatedly, and if the neighbouring subdomain is not in the 

list it will not be computed until a new list is created. It is therefore not beneficial to re 

strict the list when the number of active subdomains is low. A cut-off is therefore only cre 

ated if there are more active subdomains than an empirically chosen `MinAct` value. With 

the variable `noC` representing the total number of subdomains that the targeted platform 

can compute simultaneously, the choice `MinAct` \= max (`noC`,  max (*sxsy*,`noSched`))   
seems to give a good overall performance, and has been used in all reported results. The 

`noSched` bound assures that the list is not limited if there are much fewer active subdo-  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 13 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Algorithm**  BuildScheduleSOLAS(`L`,`CL`) 

mains than was scheduled in the previous list. The *sxsy* bound is included so that the list 

is shortened only if the current front is rather large (or thick). 

Pseudocode for creating a task schedule in the SOLAS method is given in Algorithm . 

The returned variable `noSched` is the number of scheduled tasks that have been added 

to the list `L`. 

The pseudocode for the SOLAS driver is identical to the driver SolverLAS() in 

Algorithm , except that calls to BuildScheduleLAS() are replaced by calls to the 

SOLAS-specific procedure for constructing the subdomain schedule, that is, 

BuildScheduleSOLAS() defined in Algorithm . 

**2.3 Algorithmic observations and implementation details** 

In this section we discuss implementation specific details, and some of the optimisations 

that were tested. The gain of a particular optimisation depends on its interplay with other 

imposed optimisations. In general, if an optimisation is said to give a ‘significant’ per 

formance boost, the reduction in computational time was lowered with more than %. 

From an implementation perspective, LAS is easier to implement. However, SOLAS re 

quires only small modifications of the LAS code. The workload and number of iterations 

of LAS and SOLAS differ significantly with the underlying physical problem and boundary 

conditions. A performance modelling analysis is therefore difficult to perform. 

*Initial computations* 

An expanding front passes a node only once and reaches nodes further away at later times. 

Because of this causal dependency, the parts of the solution holding larger values depend  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 14 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

on nodes with smaller values. Numerical errors introduced close to the origin will prop 

agate to the whole domain. It is therefore important that the computed solution is par 

ticularly accurate close to the initial condition. To ensure that all initiated subdomains 

fully converge, they are computed with one sweep of D PMM in the beginning of the 

algorithms. 

*Blocking for cache* 

In our initial implementation, the numerical grid was stored in row-major order. Nodes 

within one subdomain are stored in several rows, and every row of data continues into 

other subdomains. Therefore, some of the data read into cache will belong to other sub 

domains, and will not be used in computations. When data that are reused often are stored 

physically close in memory, more relevant data are loaded at once to the faster cache mem 

ory. Nodes within a subdomain were therefore ordered contiguously in memory. This op 

timisation improves the data locality, and is sometimes referred to as *blocking for cache* 

\[\]. The GPU performance was improved significantly with this optimisation, since the 

access time between global memory and shared memory differ more on GPUs than on 

CPUs. On GPUs, rows of two and three dimensional data structures are stored in memory 

allocations of ‘pitch’ bytes. The pitch is chosen by CUDA to ensue best performance when 

accessing the row addresses \[\]. In the GPU implementation, all nodes in one subdomain 

lie within one pitch so that the memory accesses are aligned for both the synchronisation 

and the computation kernels. The CPU performance was also improved, but only by a 

couple of per cent. 

*Sliding window* 

If the entire subdomain is read into the shared memory on the GPU, the performance is 

limited by the high use of shared memory. This was implemented as a first attempt, but 

led to unacceptably low speedup. Instead, we use a *sliding window* approach \[\]. When 

two adjacent nodes are updated in D PMM, they share several values in the bottom layer 

of the pyramid stencils, see Appendix B.. To reduce the number of reads from global 

memory, the entire bottom layer of nodes is first loaded to a shared memory structure, 

from where cache loads are faster. The GPU performance was improved significantly. 

A similar sliding window optimisation is used in the CPU implementation, since the 

entire bottom layer of nodes is loaded to a cached data structure. This increases the per 

formance slightly, and has been used in the numerical examples in Section . 

*List building in GPU implementations* 

Both the SOLAS and the LAS methods create lists. This is difficult to do efficiently on 

a GPU. Instead, a boolean structure indicating whether or not a subdomain should be 

present in the list is created on the GPU. This structure is transferred to the CPU where 

the list creation can be performed fast. Unfortunately, this requires the GPU to wait while 

the list is being built on the CPU. Stalling communications between the CPU and the GPU 

is known to be a bottleneck of GPU algorithms \[\]. 

*Kernel timeout on GPUs* 

In our application, we sometimes experienced difficulties launching kernels for the 

anisotropic experiments in double precision on laptop GPUs. These devices were simul 

taneously used by the operating system to display the graphical desktop environment. It  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 15 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

appeared that when launching a kernel on many subdomains, we received the error mes 

sage CUDA\_ERROR\_LAUNCH\_TIMEOUT. According to the CUDA Driver API \[\], 

this error is returned when a kernel takes too long to execute. To avoid the device making 

the operating system halt, Nvidia has implemented a specific time limit on the execution 

time of kernels. If this time limit is exceeded by the kernel, the kernel will simply not 

launch or in some cases be terminated once being launched. 

Our solution is to limit the number of launched kernels. If there are too many active sub 

domains in the list, the list is sliced into chunks on which kernels are launched. The size 

of each chunk is chosen so that the device can execute the kernel with maximum occu 

pancy. In cases where the number of active subdomains will not lead to a kernel timeout, 

the chunking is simply omitted. 

**3 Numerical verification** 

In this section, we review a series of numerical experiments conducted in order to study 

the performance of the new solvers proposed in this paper. These experiments include 

isotropic and anisotropic problems for which analytic solutions can be derived, and 

isotropic problems with changing velocity and obstacles that are impermeable or very 

slowly permeable. Our final case is that of simulating a geological fold, based on seismic 

field data. 

**Definitions** 

For any scalar field *ξ* represented by values in every grid node (*i*, *j*, *k*), we will apply the 

discrete *L* norm 

   

 *ξ* *L* \=   
*i*,*j*,*k ξ* *i*,*j*,*k*   
*N* , () 

where the summation is performed over all *N* nodes in the domain. We refer to the grid as uniform only when *dx* \= *dy* \= *dz*. 

**Computing platforms** 

The numerical experiments reported herein have been conducted on a  core CPU re ferred to as *SandyBridge*, and a Nvidia Tesla K GPU, hereafter called *K*. The capabil ities of these platforms are summarised in Tables  and , respectively. 

All computing times presented in this section are from implementations of the algo rithms, one multicore CPU variant using OpenMP \[\], and one GPU variant using CUDA \[\]. For all test cases, we report on the fastest solution extracted from at least three ex ecutions of the chosen algorithm on the chosen platform. In general, we have observed that the computing times have only insignificant variations from one run to another. 

**Sizing subdomains** 

In general, D PMM is an efficient algorithm with significant parallel capabilities. For each subdomain, one needs to store both the internal nodes and the ghost nodes to the sub domain. The domain decomposition is slightly overlapping, considering that the ghost nodes are inner nodes of other subdomains. Therefore, the smaller the subdomains get, the larger will the fraction of ghost nodes relative to inner nodes become. From a global  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 16 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Table 1 Specifications for the Intel Xeon E5-2650 2** × **8-core CPU (one of two identical NUMA** 

**domains).** 

**SandyBridge: Intel Xeon E5-2650** 

CPU frequency 2.0 GHz   
Number of cores 2×8   
Single precision 256 GFLOP/s   
Double precision 512 GFLOP/s   
Shared L3 2×20 MB   
Memory 32 GB DDR3   
Memory bandwidth 2×51.2 GB/s 

**Table 2 Specifications for the NVIDIA Tesla K20 GPU.** 

**K20: NVIDIA Tesla K20** 

CUDA cores 2496   
Single precision 3520 GFLOP/s   
Double precision 1170 GFLOP/s   
Memory 5120 MB GDDR5   
Memory bandwidth 208.0 GB/s 

perspective, a reduction of the subdomain size will increase the memory footprint and the 

communication overhead due to storing and handling a growing number of ghost nodes. 

Still, a relatively small subdomain size may offer an advantage in that the list can better 

mimic the front propagation, thus leading to a more efficient exploitation of the underly 

ing causality. Moreover, nodal values of smaller subdomains are more likely to fit in fast 

memory locations, thereby increasing the computational speed per subdomain. 

We experimented with different cubic sizes of subdomains for a range of problems. For 

the CPU implementation, the optimal subdomain size depends on the exact problem for 

mulation and the total number of nodes. A change in subdomain size have only small 

impact on the performance of the CPU implementation. We have therefore focused on 

the GPU implementation when choosing the subdomain size for our computations. Sub 

domains with  nodes lead to the best GPU performance in terms of wall clock time used   
for the solution. Since the computations include updates of ghost nodes, a total of ( \+ ) 

threads are used to compute a subdomain, see Appendix B. for details. The GPUs used 

for the analysis all have a warp size of , that is, a minimum of  identical operations 

per instruction is performed in a Single Instruction Multiple Data (SIMD) fashion \[\]. If 

the number of threads is not a multiple of the warp size, some GPU threads will be idling 

while the others performs computations. Having analysed our implementation with the 

CUDA Device Occupancy calculator \[\], we have concluded that the subdomain size of 

 allows for both the best occupancy and the shortest computing time. This subdomain 

size is kept constant throughout all experiments reported in this paper. 

**3.1 Synthetic examples** 

The LAS and SOLAS algorithms proposed in Section . have been tested on a range of 

problems that expose the solvers to several challenges. Here, we briefly summarise our 

findings by presenting five of those examples. For all these examples, we have computed 

solutions for grids with , ,  and  nodes along each axis, grouped accordingly 

in , ,  and  subdomains. That is, these grids represent from ,  to , 

,  unknowns. The computations have been conducted on the SandyBridge (CPU)  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 17 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Figure 8 Isosurfaces of the anisotropic solution with**   
**13 point sources located irregularly in the rectangular**   
**domain as boundary condition (Ex*A*).** 

Ex*A*:  point sources. 

and the K (GPU) platforms. In all experiments, the solvers stop when all subdomains 

are locked for computations, and the schedules are empty. 

The reported performance numbers bellow are based on computations in double pre 

cision. Due to specific hardware features of GPUs, a transition from CPU to GPU is far 

more attractive for computations in single precision. For the type of calculations done in 

this paper, single precision arithmetics can be used without noticeable loss of accuracy. 

This topic is further discussed in Section .. 

*.. Example with analytic solution* 

The first example uses the anisotropic formulation (), with constant velocity parameters 

*F* \= ., and *ψ***a** \= (., –., –.). As boundary condition we used  point sources 

spread out irregularly in the rectangular domain, defined by (, , ) ≤ **x** ≤ (, , ). The 

ratio between the maximum and minimum velocities, referred to as the anisotropy co 

efficient, is . for this problem. A visualization of the anisotropic solution is given in 

Figure . 

As shown for the Ex*A* case in Figure (a), the GPU implementations always outperform 

its CPU counterparts. Notice that SOLAS performs slightly better than LAS for large grids 

across both platforms. The difference between these two solvers is more prominent on the 

CPU than in the GPU-based computations.   
The efficiencya for eight and  cores is \-% for SOLAS, whereas it drops down to 

\-% for LAS. Many subdomains are active simultaneously in Ex*A*, so the additional 

ordering of SOLAS is worthwhile. These efficiency ratios are good, taken into account 

the sequential construction of schedules and that threads are synchronised between the 

ComputeSchedule() and SyncFromSchedule() procedures. Going from CPU to GPU 

architecture, we observe . and . times faster computations for the Ex*A* case when 

employing the SOLAS and LAS solver on the largest grids, respectively. We also run the 

same problem with an isotropic setting, *ψ***a** \= . For the isotropic case, LAS is instead 

slightly faster than SOLAS on large grids. Compared to the isotropic counterpart, the 

anisotropic computing time increases with a factor of \- on the CPU and \- on the 

GPU. The anisotropic stencil is significantly more complicated than the isotropic, and 

the solution dependencies in the anisotropic case increase the total number of times that 

subdomains are activated. 

In Figure (b), we have repeated the computing times from Figure (a), using logarith 

mic scales on both axes. The thin dotted lines illustrate linear relations to *N*. That is, the 

figure indicates that all the tested methods have a computational cost that scales with the  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 18 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

|  |
| :---- |

|  |
| :---- |

(a) Ex*A*: computing times. (b) Ex*A*: logarithmically scaled computing times. 

**Figure 9 Performance on the SandyBridge and K20 platforms for the isotropic test case with 13 irregularly located point sources (Ex*A*).** (**a**) Observed computing times (seconds) for 16 cores. (**b**) Observed computing times presented in logarithmic scales. 

grid size close to *O*(*N*), or even *superlinearly*, that is, the computing time seem to increas ing less than linearly with growing problem size. Similar behaviour is observed for all our numerical investigations. There are several factors that can help explain this observation. Small data sets cause more computing units to be idle, thereby reducing the observed com putational speed. When the data set grows in size, the level of parallelism increases and a better load balance between the processing units is achieved. The time needed to transfer data to and from the GPU is included in the time measurements. The overhead associ ated with (synchronous) data transfers between the CPU and GPU decreases as the data set gets larger. The cost of data transfer is therefore relatively larger for smaller data sets. Although we have not formally analysed the computational cost of the new algorithms, they will likely have linear complexity due to being semi-ordered and of sweeping type. The superlinear behaviour is thus caused by other artefacts as explained above. 

*Accuracy and convergence* This example has analytical solutions \[\], making it possi ble to study how the numerical solution converges towards the true solution. Let *hk* be the edge length of the *k*th grid in a sequence of increasingly finer grids used for numeri cal solution of the problems. Denoting the corresponding error *e*(*k*), the estimated rate of convergence is 

*p* \= log( *e*(*k*) *L* / *e*(*k*\+) *L* )   
log (*hk*/*hk*\+) . () 

For the grid sizes in this study, our convergence estimates are ., ., . and ., for both the anisotropic and isotopic results. Convergence estimates for single and double precision solutions are for all practical purposes identical. To reach first order conver gence, the amount of nodes initially assigned analytic values should be kept constant \[\]. The estimated convergence rates are therefore satisfactory. 

*.. Examples with curved characteristics* 

Iterative algorithms are sensitive to curved characteristics \[\]. It is therefore common to investigate how algorithms perform for problems with obstacles that force the character istics to turn sharply \[, , \].  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 19 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

![][image3]![][image4]![][image5]  
(a) (Ex*C*)Ex*B*: (Im)permeable walls. (b) Ex*D*: Velocity-cubes. (c) Ex*D*: Velocity-cubes cut. 

**Figure 10 Illustrations of solutions to examples with curved characterics.** (**a**) The distance field in   
presence of semi-permeable obstacles organised as nine walls (examples Ex*B* and Ex*C* ).   
(**b**) Three-dimensional visualisation of the distance field obtained by solving the chequerboard velocity   
problem (example Ex*D*). (**c**) Two-dimensional planar cut through the velocity-cubes distance field (example   
Ex*D*). 

In the test cases Ex*B* and Ex*C*, we consider a cubic volume divided by nine walls with 

small openings in either the upper left or lower right corners. The walls are one node 

layer thick, and are either semi-permeable (*F* \= .*dz* for Ex*B*) or impermeable (*F* \=  for 

Ex*C*). Outside these walls, the isotropic velocity is set to *F* \= . The fastest path from the 

starting point in a corner of the domain to the opposite corner is by zigzagging through 

the small openings in the walls. A sample solution is visualised in Figure (a) \[, \]. 

In the final example, Ex*D*, the computational domain is divided into a set of equally sized 

cubes. The velocity is constant within each cube, having the value of either  or . For a 

cube with velocity value , all adjacent cubes sharing a side with it will have a velocity equal 

to , and vice versa. We have placed in total  such velocity cubes within a computational 

domain measuring  length units along each axis. The boundary condition is one point 

source at the center of the domain. Figure (b) shows a volumetric view of the computed 

solution, while Figure (c) shows a planar cut of the solution through the center of the 

domain \[, \]. 

*Performance for the obstacle problems* The examples Ex*B* and Ex*C*, having barriers em 

bedded in the computational domains, are challenging since the front is forced to travel 

along irregularly shaped paths. Figures  and  summarise the computing times and the 

speedup factors for these two problems. 

The barriers are only one node thick for all grids, and the fraction of nodes that are obsta 

cles reduces as the grid gets finer. The front is small since it only propagates between two 

layers of walls. For larger grids, more subdomains fit between two barriers, and the par 

allel processing possibilities increases as the front activates more subdomains. Still, many 

schedules are built in the computation of these examples. The changing dynamics of the 

problems causes the compute times to change rather peculiarly with grid sizes, especially 

LAS since the extra ordering in SOLAS reduce the effects. At a large enough grid size, the 

LAS solvers changes to such an extent that the solve time decreases for a larger grid size. 

The irregular behaviour of LAS is present for both the CPU and the GPU implementa 

tions. Although the variations for each platform is at the same level when seen relative to 

the platform’s computational power, this issue is most pressing in the CPU environment 

due to generally longer computing times.  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 20 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

|  |
| :---- |

|  |
| :---- |
|  |

(a) Ex*B*: computing times. (b) Ex*B*: LAS scaling. 

**Figure 11 Performance on the SandyBridge and K20 platforms for the example with semi-permeable barriers (Ex*B*).** (**a**) Observed computing times (seconds) for 16 cores. (**b**) Observed speedup of LAS when moving from one (*t*(1)) to many (*t*(*p*), *p* \= 2, 4, 8, 16\) cores, and for 16 cores relative to the GPU. 

|  |
| :---- |

|  |
| :---- |
|  |

(a) Ex*C*: computing times. (b) Ex*C*: LAS scaling. 

**Figure 12 Performance on the SandyBridge and K20 platforms for the example with impermeable barriers (Ex*C*).** (**a**) Observed computing times (seconds) for 16 cores. (**b**) Observed speedup of SOLAS when moving from one (*t*(1)) to many (*t*(*p*), *p* \= 2, 4, 8, 16\) cores, and for 16 cores relative to the GPU. 

For Ex*B*, where the barriers are semi-permeable, LAS turns out to be the most efficient solver for large *N*, regardless of platform. For the impermeable barriers in Ex*C*, SOLAS would be the preferred method, although the behaviour comes very close to that of LAS when the grid gets large. For the  grid, the GPU versions of LAS and SOLAS compute the solution in . seconds, while the CPU implementations require about  seconds. It is noted that the semi-permeable case requires about two to three times as long computing time as the impermeable problem, regardless of solvers and platforms. 

Figure (b) shows that LAS applied to the Ex*B* problem gives a very good speedup, even on eight and  cores for which the efficiency is % and %, respectively. Moving on to the GPU, the LAS method runs . times faster than the  CPU cores. The speedup of SOLAS on eight cores is the same as for LAS, but there is a drop to % efficiency when going up to  cores. Similar to the observation we made for Ex*A*, this drop of speedup does not propagate over to the GPU, where SOLAS performs at same level as LAS. As shown in Figure (b), the computing times obtained for applying LAS and SOLAS to  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 21 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

|  |
| :---- |

|  |
| :---- |
|  |

(a) Ex*D*: computing times. (b) Ex*D*: SOLAS scaling. 

**Figure 13 Performance on the SandyBridge and K20 platforms for the example semi-permeable barriers (Ex*D*).** (**a**) Observed computing times (seconds) for 16 cores. (**b**) Observed speedup when moving from one (*t*(1)) to many (*t*(*p*), *p* \= 2, 4, 8, 16\) cores, and for 16 cores relative to the GPU. 

Ex*C* scale with the number of cores quit similar to the factors observed for Ex*B*. However, the performance gain of moving from  cores to the GPU is less for the impermeable problem in Ex*C*, where the speedup factor is in the range .-.. 

*Performance for the velocity-cubes problem* For the velocity-cubes problem Ex*D*, Fig ure  illustrates smooth behaviour of all algorithms with respect to the computing time as a function of the grid size. As for the two examples with embedded barriers, SOLAS proves to be an efficient solver. It performs consistently better than LAS, in particular on the CPU platform. On the GPU, the solution can be computed for the largest grid in . seconds, compared to the . seconds required by the  CPU cores. Both LAS and SO LAS demonstrate excellent speedup on the CPU, delivering %-% of the theoretical potential for eight and  cores. For both methods, the GPU calculations runs \- times faster than the  cores. 

*Overall performance for the synthetic examples* Combining the observations from the synthetic test cases discussed above, we conclude that the proposed methods are well suited for efficient computations on multicore CPUs and GPUs. Even for complicated problems these algorithms offer computational efficiency that can support interactive soft ware applications, also when handling finely spaced grids. This behaviour is due to the methods being designed to efficiently exploiting the causality of a propagating front. 

For anisotropic problems, the GPU implementations of LAS and SOLAS spend on the average around % of the time on computing, \-% on synchronisation, and less than % on list building. For isotropic problems, the computations take around % of the total time on the GPUs. Computations dominate even more on CPUs where \-% of the total time is spent on subdomain computations and the remaining time is used for synchroni sation. List creation is then less than .% of the total time, and is therefore negligible. As to be expected, the synchronisation of values is slightly more costly for SOLAS than for the LAS algorithm. 

Regardless of the test cases, the GPU implementations offer consistently faster compu tations than the CPU implementation running on up to  cores. The speedup factors when moving from  cores to the GPU ranges from  to  with an average of . For single  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 22 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

precision computations, this speedup factor varies from  to  with an average of . LAS 

is preferred for the isotropic or mildly anisotropic problems, whereas SOLAS is the best 

choice for problems with complicated velocity profiles. 

**3.2 Single versus double precision arithmetics** 

The difference in peak performance between single and double precision of the Sandy 

Bridge CPU assumes a perfectly vectorised code. Unfortunately, the large amount of log 

ical branching prohibits a vectorisation of the computationally dominating update step. 

By hardware design, GPUs compute substantially faster in single precision than in double 

precision. In particular, Nvidia’s Kepler architecture, such as our K platform detailed 

in Table , is able to perform three times as many operations per second in single preci 

sion than in double precision. In addition to the raw processing capabilities, values repre 

sented in single precision need only half the memory space that values stored in double 

precision need. Therefore, one can fit larger data sets into the limited memory of a GPU 

if single precision is used. Data-intensive applications will thereby get an additional im 

provement of performance due to better utilisation of fast memory and less shuffling of 

data between the GPU and the host computer. Across the range of synthetic problems 

investigated in Section ., we have observed an average speedup of . when comparing 

single and double precision computations on the K platform. In comparison, the cor 

responding speedup of . for the \-cores SandyBridge CPU is rather modest. On both 

platforms, these speedup factors are independent of the solvers. This observation agrees 

with the fact that the total time spent by the algorithms is heavily dominated by arithmetic 

operations, rather than memory transfers and synchronisation. 

Given the advantage of single precision arithmetics, it is relevant to investigate whether 

computations in single precision can provide a satisfactory level of numerical accuracy 

when solving problems of the form (). Empirically, we have observed that the single pre 

cision solution, *Tf* , is practically identical to the double precision solution, *Td*. For the 

largest grid, *N* \= , we have observed that  *Td* – *Tf*  *L* \= *O*(*p*), where *p* \= – for Ex*A*,   
*p* \= – for Ex*A* and Ex*D*, and *p* \= – for the complicated cases Ex*B* and Ex*C*. By the triangle 

inequality, we have 

 *T* – *Tf*  *L* \=  *T* – *Tf* \+ *Td* – *Td* *L* () 

≤ ( \+ *ε*) *T* – *Td* *L* , () 

where *ε* \=  *Td* – *Tf*  *L* / *T* – *Td* *L* . That is, if *ε* \< , the single precision error will be of the same order as for the computations in double precision. We have estimated *ε* empiri 

cally for Ex*A*, which have analytic solutions. Depending on the problem and the grid size,   
this entity stays in the range from – to – and is thus significantly smaller than one. 

Moreover, we have computed the ratio *β* \=  *Td* –*Tf*  *L* / *Td* *L* , which is obtainable also for   
problems that can not be solved analytically. We have then observed that *β* stays between   
– and – for Ex*A*\-Ex*D* across all grid sizes, and that *β* gets smaller as *N* increases, thus 

indicating convergence. 

Based on these empirical studies, we conclude that single precision computations pro 

vide high enough accuracy when using the algorithms proposed in this paper with first 

order stencils to simulate front propagation based on (). More subdomains are reacti 

vated for computation when using double precision arithmetics than in the case of single  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 23 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

![][image6]![][image7]  
(a) (b) 

**Figure 14 Illustrations of a fold simulation.** (**a**) A surface extracted in the CES application, with two cross   
sections through a three-dimensional offshore seismic volume. (**b**) The folded layers simulated on basis of the   
extracted surface and an isotropic formulation. 

precision. A subdomain is only reactivated when it receives a new value. Since the small 

est noticeable change in single precision is effectively –, the subdomains with smaller 

changes will not be reactivated in single precision, but only in double precision compu 

tations. The increase of reactivations contributes to the increase of time needed by the 

solvers in double precision compared to single precision. Note that a subdomain could be 

reactivated only if the new value *t*new differs from the old value *t*old more than a chosen 

threshold-limit in order to limit the amount of reactivations. In other words, to reduce the 

number of subdomain computations, a subdomain can be activated only if *t*old – *t*new \> *c* 

for some empirically chosen value *c*. In all our numerical experiments we have used *c* \= . 

**3.3 Fold simulation** 

Using the CES software, a user can easily segment structured surfaces from seismic data 

sets by autotracking. Figure (a) shows an autotracked surface in the CES software, ex 

tracted from a three-dimensional seismic volume as indicated by the transparent seismic 

cross sections. The physical domain measures , , and  meter along the *x*, 

*y* and *z* axes, respectively. Values in the nodes immediately above and below the extracted   
horizon were computed in CES for grids with  and  nodes, respectively. These val 

ues were then used as boundary conditions in our numerical experiments. Parallel folded 

layers were simulated by solving the isotropic eikonal equation for the different grids. An 

example of the computed solution is visualised in Figure (b). 

We have compared our new algorithms with the CES solver. The fold simulation solver 

in CES is based on a tracking algorithm, and is therefore sequential. When measuring 

the performance of the CES solver, we used a  core, . GHz Intel Ivy Bridge Core i- 

K CPU with MB L cache. This is the fastest platform available that can run the CES 

software. Applying the new algorithms on K and SandyBridge, we have measured the 

computing times for both single and double precision arithmetics. Because of the sizing 

of subdomains, the grids used by the new algorithms are slightly larger than the grids used   
by the CES solver. The new grids consist of  and  nodes, respectively. These grids 

represent from , to more than  million unknowns. Table  shows the computing 

time for double and single precision arithmetics, together with the speedup achieved by 

replacing the Ivy-based CES solver. The speedup factors are prefixed with the symbol ‘×.’  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 24 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Table 3 Comparison of computing times (in seconds) for the CES solver on Ivy and LAS and** 

**SOLAS solver on SandyBridge and K20, including the speedup factors (‘**×**’) relative to CES.** 

**Solver Device Double precision Single precision**   
 **  ** 

CES Ivy 27.68 ×1.0 354.9 ×1.0 24.91 ×1.0 319.42 ×1.0   
LAS K20 0.87 ×32.0 6.2 ×57.3 0.44 ×57.2 2.99 ×107.0   
SOLAS K20 0.84 ×33.1 6.0 ×59.6 0.43 ×58.6 2.98 ×107.2   
LAS SandyBridge 5.22 ×5.3 36.0 ×9.9 4.31 ×5.8 28.46 ×11.2   
SOLAS SandyBridge 4.66 ×5.9 32.8 ×10.8 3.81 ×6.5 26.67 ×12.0 

The CES time for single precision computations has been empirically estimated as 90% of the CES time for double precision. 

Reviewing the results, we observe that the difference between solution times for LAS 

and SOLAS are negligible when run on the same platform and with the same arithmetic 

resolution. We also observe that the speedup is higher for the GPU platforms than for the 

CPU platforms. As noted in Section ., there is no difference in solution quality between 

single and double precision computations. However, the computing times are consider 

ably lower in single precision computations, especially on the GPU. The computing time 

for the GPU is roughly halved when moving from double to single precision computa 

tions. For the largest grid in single precision, the K implementation of SOLAS is more   
than  times faster than the CES solver on Ivy.b We do not expect the speedup to be as 

impressive in all fold simulations. However, the user experience is significantly improved 

when fold simulations are performed with the new algorithms. Even on the largest grids, 

the computational time is a few seconds instead of minutes when LAS or SOLAS uses 

K instead of the current solution method. Based on this and several other numerical 

experiments with seismic field data, the proposed algorithms are now being implemented 

in the industrial code. 

The task of simulating a folded volume is rather different from the synthetic examples 

discussed initially. For fold simulations, the initial geological surface defining the bound 

ary condition will contain a significant number of subdomains. Because of this, the list 

used for subdomain scheduling will be quite long already from the start. When this list is 

long, it is easier to make full use of the available computing resources, thus resulting in a 

substantial effect of parallelisation. For simpler problems, where the boundary condition 

is represented by one or a few point sources, only a few subdomains will be scheduled as 

active. Thus, some of the available computing resources will be idling until later stages of 

the solution process when the lists grow longer. Since GPUs have an intrinsic capability 

of processing a large number of subdomains simultaneously, the simulation of geological 

folds based on LAS or SOLAS is particularly suitable for GPUs. 

**4 Conclusions** 

Two new parallel solvers are presented, together with several numerical experiments con 

ducted on both multicore CPUs and GPUs. All solvers use a domain decomposition ap 

proach, where each subdomain is surrounded by a layer of ghost nodes. A subdomain 

needs to be processed only if its local boundary condition has changed since the previ 

ous computation. Due to this observation, a computed subdomain is locked for further 

computations until its boundary condition changes, and the subdomain is reactivated. 

The presented algorithms differ mainly in the way subdomains are scheduled for com 

putations. LAS places all active subdomains in a list, and updates the subdomains in this  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 25 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

list in parallel. The GPU implementation of LAS builds the lists on the CPU, but performs 

all computations on the GPU. SOLAS is very similar to LAS, but creates the lists slightly 

differently to enforce a more causal ordering of the subdomain computations. 

From empirical studies, the sequential solution times seem to scale linearly with the 

total number of nodes. The parallel solvers appears to scale somewhat superlinearly, ac 

credited to the increase in parallel processing efficiency obtained for larger grids. Some 

algorithmic details are crucial for achieving good performance. For instance, if a schedule 

is only used once by LAS and SOLAS, many more and longer lists are built. For anisotropic 

cases, repeated activations between neighbouring subdomains are more common than for 

isotropic cases. It is therefore important not to enforce a too strict ordering of subdomains. 

For some examples the LAS ordering is not very different from the SOLAS ordering, and 

the extra overhead of the SOLAS ordering is then not worthwhile. 

We have investigated the potential for accelerating an industrial software for simula 

tion of geological folding by introducing new algorithms for numerical solution of static 

Hamilton-Jacobi equations. The computing times are reduced from several minutes to 

seconds, enabling the software to be used interactively even for large three-dimensional 

data sets. 

**Possible extensions of the presented work** 

In our implementation, the computing and synchronisation procedures are entirely com 

pleted before the list building starts. Future implementations can overlap the building of 

new schedules with the processing of previous schedules. That is, one processor can build 

the new list while other processors reuse the old list. This will result in reduced overhead 

of list building and processor synchronisation. On the GPU, the copying of data to and 

from the device can be done asynchronously with the computations, and the utilisation of 

the CPU will increase. 

The ghost nodes assure that subdomains can be updated simultaneously by different 

processors. No data is communicated between subdomains during the update procedure. 

Therefore, a parallel implementation iterates through the grid in the same way as a sequen 

tial implementation. The computed solutions are independent of how many processors are 

used, and the computed solution is always of good quality. 

Ghost nodes carry the drawback of needing a synchronisation procedure. Before a 

synchronisation of nodal values, all computations must be finalised and the processing 

threads must be synchronised. The pHCM method avoids synchronisations of threads by 

allowing memory to be shared between different threads \[\], resulting in so called *false* 

*sharing*. Because of the stable upwind stencils, the pHCM method converge to the cor 

rect solution. However, the algorithmic behaviour changes noticeably depending on by 

which processor, and in which order, memory is accessed. Similarly to pHCM and FIM, 

the ghost nodes can be avoided in LAS and SOLAS. Some minor changes of the reacti 

vation system is then needed, and the synchronization can be entirely avoided. Values are 

instead ‘synchronised’ between adjacent subdomains directly when computed. We expect 

the performance when running on a few threads to improve significantly, but the number 

of iterations will change with the number of threads. 

Other methods than D PMM should be investigated for the nodal updates within each 

subdomain. It is difficult to foresee which type of methods would be better suited. Itera 

tive methods such as SOFI, FIM and Locked Sweeping \[\] might perform well since the  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 26 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

![][image8]  
(a) Stencil tetrahedron. (b) One stencil shape. (c) Independent updates. 

**Figure 15 Illustration of numerical stencil and nodal update procedure of 3D PMM.** (**a**) A stencil shaped   
as a tetrahedron, with four characteristics marked with *A*, *B*, *C* and *D*. Only the solution along characteristic *D*   
is acceptable. (**b**) The pyramid-shaped stencil used in 3D PMM, consisting of eight tetrahedrons. Nodal values   
from the bottom layer is used to update the value of node *Ti*,*j*,*k* . (**c**) A layer of pyramid stencils. All top nodes in   
the layer can be updated in parallel, since none of the computed values affects the values of other nodes   
being updated. 

subdomains are rather small. Second order stencils can be implemented with almost no 

algorithmic alterations, thereby increasing accuracy further. 

**Appendix 1: Efficient formulation of the eikonal ‘pyramid’ stencil** 

The most common way to propagate a front numerically is through use of semi 

Lagrangian \[\] or conditional upwind \[, \] discretisations. The semi-Lagrangian ap 

proach searches through all possible characteristics for the one yielding the minimal dis 

tance. The conditional upwind stencils first solve for a solution, and thereafter accepts the 

new estimate only if the associated characteristic originates from within the spatial grid 

element, defined as the convex hull of the nodes supporting the stencil. Both approaches 

are identical on eikonal formulations but may differ in anisotropic cases \[\]. If the char 

acteristic curve originate from outside of the element, the new estimate is not based on 

upwind values, and should not be accepted. Due to the isotropic nature of the eikonal 

problem one can often formulate upwind conditions for the acceptance of a new estimate 

\[\]. These conditions can significantly reduce the amount of computations needed to 

update a nodal value. 

In this appendix we present an efficient stencil for the eikonal equation, derived as a 

conditional upwind stencil. The stencil element is shaped as a tetrahedon defined by the 

corner nodes *Xi*,*j*,*k* with values *Ti*,*j*,*k* \= *T*(*Xi*,*j*,*k*), as visualised in Figure (a). The value *T*,, 

for the top node is to be updated using the values *T*,,, *T*,, and *T*,,. Four character 

istic curves (dashed arrows) are shown with their entrance points on the *bottom* surface 

defined by nodes *Xi*,*j*,. Only the solution associated with the characteristic curve *D* would 

be acceptable since the other curves originate from outside the stencil element. When up 

dating a node in D PMM, the tetrahedron shape is used eight times in a configuration 

that resembles a *pyramid*, shown in Figure (b). 

In line with general consensus, the velocity *F* is assumed to be constant within the ele 

ment. When updating *T*,,, as in Figure (a), only nodes with strictly smaller *T* values 

than *T*,, are used. A larger value implies that the front moves away from *X*,,, and to 

wards the position of the larger value. Initially we assume that no new value is found, that  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 27 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

is 

*t*new ← *T*,,. () 

From the characteristic equations () we know that the characteristic coincides with the 

gradient **n**, and since **n** \= ∇*T*/ ∇*T*  \= *F*∇*T*, we get the following estimates 

**n***x* ← *FT*,, – *T*,,   
*dx* , **n***y* ← *FT*,, – *T*,,   
*dy* . () 

These estimates and the assumption of constant velocity correspond to modelling an ar 

riving front as planar \[\]. If **n***y* \> , the ray reaching *X*,, cuts the bottom surface outside 

of the tetrahedron, in areas with negative *i* offset (ray *A*). The shortest distance from within 

the tetrahedron then originates from the side with *i* \= , that is, the new estimate is given 

from a two-dimensional diagonal stencil solution `s`*D*(*X*,,,*X*,,) \[\]. A similar obser 

vation holds for **n***x*; in order for the characteristic to cut below the line connecting *X*,, 

and *X*,,, the condition **n***x dy* \< **n***y dx* must hold. Otherwise (ray *C*), the minimal new es 

timate is `s`*D*(*X*,,,*X*,,). These conditions results in the following upwind conditions 

`if` (**n***y* \> ) *t*new ← min  *t*new,`s`*D*(*X*,,,*X*,,) , ()   
`else if` (**n***x dy* \> **n***y dx*) *t*new ← min *t*new,`s`*D*(*X*,,,*X*,,) . () 

Next, **n***z* is estimated using the unit argument of the normal, **n***z* \=– **n***x* – **n***y* . To make sure   
that the ray originates from a location with *i* ≤  (to remove ray *B*), the following upwind 

condition is needed 

`else if`  **n***z dx* \< **n***x dz*  *t*new ← min  *t*new,`s`*D*(*X*,,,*X*,,) . () 

This condition also assures that **n***z* \> . If all upwind conditions are true, the characteristic   
curve that cuts node *X*,, enters the tetrahedron through its base (ray *D*), and the new 

arrival time estimate is 

 

`else`   
*t*new ← *T*,, \+ **n***zdzF* . () 

**Details for the anisotropic stencil** 

The anisotropic stencil for equation () uses a direct conditional upwind discretisation in our implementation. Note that nodes upwind of the new solution point can be used in the anisotropic stencil, as long as the point where the characteristic curve enters the stencil lies upwind of the updated node \[\]. 

**Appendix 2: Subdomains** 

Subdomains are surrounded by a layer of virtual nodes referred to as *ghost nodes*. These ghost nodes are copies of nodes in adjacent subdomains. Since each ghost node exists at multiple locations, their values must be kept equal at all locations. We refer to the pro cedure of comparing such nodal values as synchronisation. Figure  is a D illustration  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 28 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Figure 16 Two neighbouring subdomains, without**   
**and with ghost nodes, top and lower illustrations.**   
Ghost nodes are coloured in green. Columns of nodes   
that exist in both subdomains have been marked with   
waves or stripes. Arrows show the direction in which   
nodal values are copied between the subdomains.   
Notice that two layers of nodes are synchronised   
between subdomains. 

of subdomains in their logical grid (top) and with the ghost nodes marked by light green 

(bottom). Columns of nodes that exists in both subdomains have been marked with waves 

and stripes in the lower figure. 

Nodes belonging to a subdomain, including the ghost nodes, will during computations 

only be accessed, and written to, by the thread(s) updating that subdomain. Two subdo 

mains can be computed simultaneously without any risk of memory interference since all 

data is private to the updating thread(s). Ghost nodes have been used and discussed in 

relation to front propagation problems in other works \[, \]. Our use is different, and 

details regarding the synchronisation of ghost nodes are given in Section B.. The compu 

tation of a subdomain using D PMM is discussed in Section B.. 

**B.1 Synchronisation of nodal values** 

Consider the two-dimensional case when comparing values from nodes in subdomain 

(*ii*, *jj*) to values in subdomain (*ii* \+ , *jj*), as shown with the two lower arrows in Figure . 

When a computed nodal value, *ts* in (*ii*, *jj*) is compared to its copy *tt* in (*ii* \+ , *jj*), the copy 

is only changed if its value is larger than the computed nodes, that is *tt* \= min (*tt*,*tc*). If 

at least one node of the inner part of (*ii* \+ , *jj*) (lower solid arrow in Figure ) receives a 

new value, the subdomain is activated by releasing the compute lock, `CL`(*ii*\+, *jj*) ← `Open`. 

A subdomain that receives new values only on its out-most boundary (dashed arrow) does 

not need to be activated. 

In SOLAS we need the minimum nodal value that activates a subdomain. The minimum 

activation value is computed efficiently in our GPU implementation using parallel reduc 

tion. Note that if a tetrahedron (diagonal) stencil is used, the corners will be shared and 

updated by eight (four) blocks in three (two) spatial dimensions. However, by synchronis 

ing values in one direction at a time, the corner nodes receive correct solution values. 

**B.2 Computation of a subdomain** 

In all our numerical experiments, D PMM has been used for the computation of new 

nodal values in the subdomains. Here we give a brief description of the method \[, \]. 

A specific *pyramid* shape of the stencil is used, as shown in Figure (b), where the 

value *Ti*,*j*,*k* for the top node of the pyramid is updated. D PMM computes *layers* of nodes, 

as illustrated in Figure (c). Nodes in a layer have one index (*i*, *j* or *k*) in common. All 

nodes in one layer can be updated entirely in parallel since there are no dependencies 

between the computed nodes, since the value of a top node does not depend on other top  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 29 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

**Algorithm**  *i*SubSweep(`sd`,`dir`) 

**Algorithm**  ComputeSubdomain(`sd`) 

nodes, see Figure (c). After processing a layer, the computed values are used as bottom 

pyramid values to compute the next layer of nodes. An iteration through the domain for 

an increasing or decreasing common index is referred to as a *subsweep*. A set of all six 

subsweeps, increasing and decreasing *i*, *j*, *k* indices, are referred to as a *sweep*. 

Conceptually, the ghost nodes in a subdomain define the boundary condition of the sub 

domain. In our implementation we compute new values also for the ghost nodes. When 

the ghost nodes are updated, no values from other subdomains are used in the computa 

tions. Therefore, ‘half’-pyramids are used at the sides, and ‘quarter’-pyramids in the cor 

ners. Pseudocode for a subsweep in the *i* direction is given in Algorithm , where the 

layers defined by *i* \=  and *i* \= *bx* \+  contains only ghost nodes. 

A subdomain that is activated receives (at least) one surface with two layers of nodes in 

the synchronisation. Assuming that these ‘thick’ boundary values are correct, the charac 

teristic curves are extrapolated to the rest of the subdomain in the sweeps as new solution 

values are computed. These ‘thick’ boundary conditions explain why the solution con 

verges fast. On a local scale, the characteristic curves are straight lines, or only slightly 

curved, and the solution within a subdomain will often converge in only one sweep. The 

final solution remains the same in all our examples whether one, two, or more sweeps are  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 30 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

used when computing subdomains. Since the total computational time is minimal when 

only one sweep is used, we have used only one sweep in the numerical experiments. 

After processing a subdomain, the compute lock is closed. As explained above, it can be 

opened again during a synchronisation. Algorithm  presents pseudocode for the com 

putation of new nodal values in subdomain `sd`, using D PMM. This algorithm, named 

ComputeSubdomain(), is called from the procedure in Algorithm . 

**Competing interests**   
The authors declare that they have no competing interests. 

**Authors’ contributions**   
The main idea of this paper was proposed by TG, based on earlier research by ØH. TG and MS contributed with the   
implementation and numerical experiments. TG and AM performed the analysis, and prepared the manuscript initially. All   
authors contributed to the manuscript, and have read and approved the final manuscript. 

**Acknowledgements**   
The presented work was funded by Statoil ASA through the Akademia program, and the Research Council of Norway   
under grant 202101/I40. The work has been conducted at Kalkulo AS, a subsidiary of Simula Research Laboratory. 

**Endnotes** a Here, the term ‘efficiency’ refers to the ratio between observed speedup and the theoretically possible speedup. For   
instance, the LAS speedup of 6.0 for Ex*A* running eight cores gives an efficiency of 6.0/8 ≈ 75%.   
b The CES software is only available on the Ivy platform running in double precision. Experiments indicate that single   
precision computations on Ivy runs at about 90% of the time used in double precision. Therefore, the CES times for   
single precision experiments are constructed by multiplying the observations for double precision with a factor 0.9. 

Received: 24 February 2014 Accepted: 30 June 2014 Published: 24 Jul 2014 

**References**   
1\. Sethian JA: *Level Set Methods and Fast Marching Methods: Evolving Interfaces in Computational Geometry, Fluid*   
*Mechanics, Computer Vision, and Materials Science*. 2nd edition. Cambridge: Cambridge University Press; 1999\.   
2\. Rawlinson N, Hauser J, Sambridge M: Seismic ray tracing and wavefront tracking in laterally heterogeneous   
media. *Adv. Geophys.* 2009, 49:203-267.   
3\. Gillberg T: A semi-ordered fast iterative method (SOFI) for monotone front propagation in simulations of   
geological folding. In *MODSIM2011, 19th International Congress on Modelling and Simulation*; 2011:641-647.   
4\. Gillberg T: Fast and accurate front propagation for simulation of geological folds. *PhD thesis*. University of Oslo;   
2013\.   
5\. Natvig JR, Lie K-A: Fast computation of multiphase flow in porous media by implicit discontinuous Galerkin   
schemes with optimal ordering of elements. *J. Comput. Phys.* 2008, 227(24):10108-10124.   
6\. Kornhauser ET: Ray theory for moving fluids. *J. Acoust. Soc. Am.* 1953, 25(5):945-949.   
7\. Podvin P, Lecomte I: Finite difference computation of traveltimes in very contrasted velocity models: a massively   
parallel approach and its associated tools. *Geophys. J. Int.* 1991, 105:271-284.   
8\. Li S, Mueller K, Jackowski M, Dione D, Staib L: Physical-space refraction-corrected transmission ultrasound   
computed tomography made computationally practical. In *Medical Image Computing and Computer-Assisted*   
*Intervention \- MICCAI 2008*. Berlin: Springer; 2008:280-288. \[*Lecture Notes in Computer Science*, vol. 5242.\]   
9\. Wallman M, Smith NP, Rodriguez B: A comparative study of graph-based, eikonal, and monodomain simulations   
for the estimation of cardiac activation times. *IEEE Trans. Biomed. Eng.* 2012, 59(6):1739-1748.   
10\. Jeong W-K, Fletcher PT, Tao R, Whitaker R: Interactive visualization of volumetric white matter connectivity in   
DT-MRI using a parallel-hardware Hamilton-Jacobi solver. *IEEE Trans. Vis. Comput. Graph.* 2007, 13(6):1480-1487.   
11\. Rawlinson N, Sambridge M: Multiple reflection and transmission phases in complex layered media using a   
multistage fast marching method. *Geophysics* 2004, 69(5):1338-1350.   
12\. Petersen SA, Hjelle Ø, Hustoft S, Haubiers M: Process based data-restoration and model-reconstruction workflow   
for seismic interpretation and model building. In *EAGE 74th Conference & Exhibition, Extended Abstracts*; 2012\.   
13\. Petersen SA, Hjelle Ø: Earth recursion, an important component in shared earth model builders. In *EAGE 70th*   
*Conference & Exhibition, Extended Abstracts*; 2008\.   
14\. Cerveny V: *Seismic Ray Theory*. Cambridge: Cambridge University Press; 2001\.   
15\. Hjelle Ø, Petersen SA: A Hamilton-Jacobi framework for modeling folds in structural geology. *Math. Geosci.* 2011,   
43(7):741-761.   
16\. Gillberg T, Hjelle Ø, Bruaset AM: Accuracy and efficiency of stencils for the eikonal equation in earth modelling.   
*Comput. Geosci.* 2012, 16(4):933-952.   
17\. Chopp D: Recent advances in the level set method. In *Handbook of Biomedical Image Analysis*. Edited by   
Micheli-Tzanakou E, Suri JS, Wilson DL, Laxminarayan S. New York: Springer; 2005:201-256.   
18\. Tsitsiklis JN: Efficient algorithms for globally optimal trajectories. *IEEE Trans. Autom. Control* 1995, 40(9):1528-1538.   
19\. Cacace S, Cristiani E, Falcone M: Requiem for local single-pass methods solving stationary Hamilton-Jacobi   
equations? 2013 \[arXiv:1301.6775\]   
20\. Sethian JA, Vladimirsky A: Ordered upwind methods for static Hamilton-Jacobi equation: theory and algorithms.   
*SIAM J. Numer. Anal.* 2003, 41(1):325-363.  
Gillberg et al. *Journal of Mathematics in Industry* 2014, 4:10 Page 31 of 31 http://www.mathematicsinindustry.com/content/4/1/10 

21\. Hjelle Ø, Petersen SA, Bruaset AM: A numerical framework for modeling folds in structural geology. *Math. Geosci.*   
2013, 45(3):255-276.   
22\. Zhang J, Huang Y, Song L-P, Liu Q-H: Fast and accurate 3-D ray tracing using bilinear traveltime interpolation and   
the wave front group marching. *Geophys. J. Int.* 2011, 184(3):1327-1340.   
23\. Herrmann M: A domain decomposition parallelization of the fast marching method. In *Annual Research Briefs 2003*.   
Stanford: Center for Turbulence Research, Stanford University; 2003:213-225.   
24\. Tugurlan MC: Fast marching methods \- parallel implementation and analysis. *PhD thesis*. Louisiana State University   
and Agricultural and Mechanical College; 2008\.   
25\. Zhao H-K: A fast sweeping method for eikonal equations. *Math. Comput.* 2004, 74(250):603-627.   
26\. Qian J, Zhang Y-T, Zhao H-K: A fast sweeping method for static convex Hamilton-Jacobi equations. *J. Sci. Comput.*   
2007, 31(1-2):237-271.   
27\. Zhao H-K: Parallel implementations of the fast sweeping method. *J. Comput. Math.* 2007, 25(4):421-429.   
28\. Detrixhe M, Gibou F, Min C: A parallel fast sweeping method for the eikonal equation. *J. Comput. Phys.* 2013,   
237:46-55.   
29\. Weber O, Devir YS, Bronstein AM, Bronstein MM, Kimmel R: Parallel algorithms for approximation of distance maps   
on parametric surfaces. *ACM Trans. Graph. (TOG)* 2008, 27(4):104:1-104:16.   
30\. Gillberg T, Sourouri M, Cai X: A new parallel 3D front propagation algorithm for fast simulation of geological folds.   
*Proc. Comput. Sci.* 2012, 9:947-955. Proceedings of the International Conference on Computational Science, ICCS   
2012\.   
31\. Gillberg T, Hjelle Ø, Bruaset AM: A parallel 3D front propagation algorithm for simulation of geological folding on   
GPUs. In *EAGE 74th Conference & Exhibition, Extended Abstracts*; 2012\.   
32\. Gremaud PA, Kuster CM: Computational study of fast methods for the eikonal equation. *SIAM J. Sci. Comput.* 2006,   
27(6):1803-1816.   
33\. Hysing S-R, Turek S: The eikonal equation: numerical efficiency vs. algorithmic complexity on quadrilateral grids.   
In *Proceedings of ALGORITMY*; 2005\.   
34\. Chacon A, Vladimirsky A: A parallel heap-cell method for eikonal equations; 2013 \[arXiv:1306.4743\]   
35\. Chacon A, Vladimirsky A: Fast two-scale methods for eikonal equations. *SIAM J. Sci. Comput.* 2012, 34(2):A547-A578.   
36\. Glover F, Klingman D, Phillips N: A new polynomially bounded shortest path algorithm. *Oper. Res.* 1985, 33(1):65-73.   
37\. Dejnozkova E, Dokladal P: A parallel algorithm for solving the eikonal equation. In *IEEE International Conference on*   
*Acoustics, Speech, and Signal Processing, 2003*. *Volume* 3; 2003:325-328. IEEE.   
38\. Jeong W-K, Whitaker RT: A fast iterative method for eikonal equations. *SIAM J. Sci. Comput.* 2008, 30(5):2512-2534.   
39\. Fu Z, Jeong W-K, Pan Y, Kirby RM, Whitaker RT: A fast iterative method for solving the eikonal equation on   
triangulated surfaces. *SIAM J. Sci. Comput.* 2011, 33(5):2468-2488.   
40\. Bak S, McLaughlin J, Renzi D: Some improvements for the fast sweeping method. *SIAM J. Sci. Comput.* 2010,   
32(5):2853-2874.   
41\. Jeong W-K, Whitaker RT: A fast iterative method for a class of Hamilton-Jacobi equations on parallel systems.   
*Technical report*. University of Utah; 2007\.   
42\. Tanenbaum AS: *2\. Modern Operating Systems*. 3rd edition. Upper Saddle River: Prentice Hall; 2007\.   
43\. Kjolstad FB, Snir M: Ghost cell pattern. In *Proceedings of the 2010 Workshop on Parallel Programming Patterns*. ParaPLoP   
’10. New York: ACM; 2010:4-149.   
44\. Rivera G, Tseng C-W: Tiling optimizations for 3D scientific computations. In *Proceedings of the 2000 ACM/IEEE*   
*Conference on Supercomputing*. Supercomputing ’00. Washington: IEEE Computer Society; 2000:1-23.   
45\. Nvidia: CUDA C programming guide; 2012 \[http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html\]   
46\. Michéa D, Komatitsch D: Accelerating a three-dimensional finite-difference wave propagation code using GPU   
graphics cards. *Geophys. J. Int.* 2010, 182(1):389-402.   
47\. Nvidia: CUDA C Best Practices Guide; 2012 \[http://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html\]   
48\. Nvidia: CUDA driver API; 2013 \[http://docs.nvidia.com/cuda/cuda-driver-api/index.html\]   
49\. OpenMP.org: OpenMP application program interface; 2013 \[http://openmp.org/\]   
50\. Nvidia: What is CUDA; 2013 \[https://developer.nvidia.com/what-cuda\]   
51\. Nvidia: CUDA occupancy calculator; 2012\.   
52\. Kadlec B, Dorn G: Leveraging graphics processing units (GPUs) for real-time seismic interpretation. *Lead. Edge*   
2010, 29(1):60-66.   
53\. Cristiani E, Falcone M: Fast semi-Lagrangian schemes for the eikonal equation and applications. *SIAM J. Numer.*   
*Anal.* 2007, 45(5):1979-2011.   
54\. Sourouri M: A parallel front propagation method: simulating geological folds on parallel architectures. *Master’s*   
*thesis*. University of Oslo; 2012\. 

10.1186/2190-5983-4-10   
**Cite this article as:** Gillberg et al.: Parallel solutions of static Hamilton-Jacobi equations for simulations of geological   
folds. *Journal of Mathematics in Industry* 2014, 4:10

[image1]: <data:image/png;base64,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>

[image2]: <data:image/png;base64,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>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIUAAACFCAYAAAB12js8AAB4bklEQVR4Xuzd9YNV1fc//oOiSMfQrSKdUirdDN3dpTRId3cLE0zfuTU359650zB0N4iK3Qq2oigxs77PtdaE78/L73/ADw9m5nLj3L3X2Xvtc/bZx9i5Zg2tDw4WS3N/7skVAZEQMyDXiGByTw2m8Mkqc14wZSwMpi0TVfrSYDqyAo+tVJmr8PzFwWSdqcLHBdOKIcFkn6RiZ+D95wZT6MJceK15K14bpY7GBVMgJpgOblUZCXhPJ7hVhieY0vB35A4VtyWY3lmC12xSsRuD6fDmYNq7UZki8L32BNOsWWrvXmxzGt4rvQD/HXIomN56S+3k9w4NprDdKgp/W/E5nvUqFt8zDp9h2qtcpmBKxXZ541WmC7DdEXuUKwx/4ztlRqnYbXh/lFPEYhW7AO+Dz/VPVuuH4rEJeN0sFYMyz1ytZStQZunLgmndBJX2drDUS8xU5ZqCbR6FbRyoYlGfh2F7ro2wIljrnR1CPBiL8MtSwxCTYDmsyxUOzkIGuUorXy2D3E0M+r2vujfAoPgOBn02Sj1+06AfZ+A1/ZQFz8kahue0V+bmBj1oa5C1oUrGa9MGGpQ4QTnn4Odqg54E1MUIgxzbDfrrrHp83aDbaQYKX3nj8PwovNcWFbsI7/e2QfuGqrSl+Ay83rJbOc3Y/kSDfnuiLuD9fEnYtjhlMeE9fQadumjQL49UYrJBngS8V6xKwfu5Fhp0YITyzMbnzsN7b1Up0XgPvE/YVnUZ7/f4hkH/XFA+bPfJQ3gsWaVvxjYtMMg/UaUOwnM6GpTQWD1CecWjzO3tFZenvb9Btj7qu2l4Dsr96hBl6WTQFyj7B/2Voym25yU8t5wyP2PQIdTrjlxrYBosy7UJ8fA0KJ4Gxf8GxTL8w8HANsFGOJwrFsyFDYqroOLroEBaIhj6K1tnfNHpBn0/TcUjSI6MM+jPlSp7AyoFfwf6qxQEQVQLg451UDH48nYOiikqsMSgS+F4n43qXSfe45pB72Uoc5hBt04alPNjrnsGncHj/kiVscegkKl43gIV2ITg4YAJUV43Cv0ovpdLXfvUoLuPDfrpofrhAb4XguXdr/B93Sr1OAoWQRF1UJl2GHR0V8FnhOPz0hBw7nCV5MB2/VCwjR8hEOyHEeA+xd/nE/y0bFZXEPg+fG/XVBWP8ojvatDZzioG5ZWOckrsqxKGa7n+vVqdmKjl/vkU9Qg7pa+HQR/3VYHX8F3rYlsrKfvzCF7Ua3SuzbABpuTirsTgfmQF/5EbEFsgPlciuPEm9sq5GmDvwgbbOqm/EAhXudJHqodrDLq/DoU0Q5nGG5SJKLb0V25sbAJen9RHncNr3GNQ+NOVBYXjRqHnXFFPsCdnWFA4qYqD4BdUmNuhbHZsIyrOHqHM2IvTsCfGrle2bQjEGHyWR9nxmji/QXceKv8lrfj0M+oEPjMGAeA9Z9D7D5QFLZPLi+DwqyS0TDZ8TuwGlYqAiMZPV4RKTdLPsVvVvc80SK4dV95obfFyLisvAsyyHO85WyWhkk+MKCgjRzf8fxf8Xz+ViXKyTirYkX5di0CAzHHqFNcDyj6hs3q/N4IArXJiVeV6AUGDejXlWpVb5/yTcY5hcEK5Hn+wiNwnHsgVAFdxbFwN5WqEL41g+HCk+haR7cfPx2vUSXwpD4Lg+50qG4WUgr0pdYoyIXIj8MWco1QABRCYpM0886G5PRGH155TJuyZd84jOH5QWen4PDS5X/2ifiFUJgra71ZJKPAUNPHOwypiHwoCFWSOVzYEhvs0KpNbCzjxLVqgbIOO31GZt1CI+Ps0gi/+uHKi2zK5ECwxKhV7eRTe1x6uAnjsqFm7GJGk2/XdfZWOz0lBMD26p75CQMThe311St30aRcXWKzcKJNUlIl3tOLyMvfTx5j7LS3XH3erpFlonfHYo3XqKALjIzzvu3EqCfXlb4ZtqK3cJfG+qNf9uSwQaWgvwfY+DYqnQfGfQcHDzrwcIs7QpiU5lw1Jprs03vhFlfIqfqL5f/KWMqNL+GsVCne2urkaG4bKuLlVufiLLkOhTlDc9MVPxmdNUT4ETBqeY12ofCicKyjAgFU9/lq7jKNp6voHaI5RaefeV2YkjS5UQoJdnUUlJqGpD/hULAIkDBWXcURFoNIO4PlZf6poVMpevO7KX+oGupQQvGd4Jrb1d3UQgRiDYAqkqzB8vxi8rz9R8ecdz8DnW1VSmgbRieuKc5QP8T0SPYq7kiffIFhc6gK204Pv7VymMpZouUROVbEorzTuXkerjJX4HJTZpc3qIbbnE/xMma4eoiuxIMd4xIECftSXvzUee0klBOFzntUdnnF9c2DwoILxkNXg4xA8wmCcVHIOwa0Dc5dB4VVBRb6k3HhzR1fNeNkZRPVvaCES31SPsYEOjABuHVRPUAEZ6NdPrFacZbtW4L2WKTs/hmTQvVm58BovCjzn21wIiB+/1NaB/ZiD16Fg066qD1DgNlRg+hFlQv6RgBbD5lUmp0GhaCFinWq/DX+fMGh3pkrhiv8Vv5sVB0TKXWwzAmTfMRVyUgMnPlGF4P3iOJ9JVHZsTzw+N+ukiuGcBduZ9b6Ktmoye+y0+gyBnf1jwXdMxjZ5w1BuG5UbElCmpuXKifJyY0R1Zq3yQXYArcFhZcbO9A9+ps9T3yForiE43p+gvD3wHNRb5MvqeDW8Z9l/1XER/P9z+OxCysJBwQeluPJZXEVsEILAXzMXWgcTgiG9sYpticd64c36q28RAKfQfH2zXXHrcPsQAsWjUvZgo/B4Ui4fAsSLJi92m7IeQKEiEFKilBm/X0LC9/l76hG6CAsK/Nvf1TkUaDqy91uPlAktgx3Pj3IrDoAIVH50snoHe+IWE5rFSMUtwEY8bvtZrcDzk/EZMWfVflRo5t8YluM18d+rLXhNCFqIzVG54vC8hILPiMDnhuF9wu3KeUa7pmt/q1OfoJW7hJHNQ2XCc3/B+373ibp2Ad8RQZEcpRw8UkK5xGxX7l0oux1oTbaoZHTLPpThDy71FXaik+g2ftmnkhEQfyAwzH3VEdRTCOot0ERxfSZyy19TOVDfpooaKMw64GlQPA2K/wqK6BHB5KttCHMdQ4adnFAyziG4y4htpdI640PQXyUMVE/W403QhTwOVzZ0Dw/RbCduV3+eR+Xh7/h9yoPnJKFQbDEqJgIbhYJ05UGFJqcZMgRkAVTS1ZsG/ZSj7Pj7XRR0/Enl/xTbhc8IC6gtqLDdKKhtHrUdw9jVeM8FYWoJCnA1nnfgY2VFYBxAd3IyW52GMFRgzOf4+Y1ag+cvR9eyIFStREBsTyn4DP48DpRQv4rEdmXcxc909d4TBKy94Dv8gOAORcWbrSoV3Z4L3ZErQSWgO4lDwmzJFcBnu7Ht1nfUJTznIbYx72AZd9GR2DGzI1U86uMJuqCU4eoyZHQqqEMX6jMVQeJqrHi4GvsK6v9FZUE8GHwuw4MIYq5WehwivbMKdNUcwt9bcUD8heQna6hKG4/nTMGHvq18aDkysKHpu9V1VIgfFZHuVmdRYImobLtLpWFPP30RrUeK4qTxOL6ww6s4m3dCCvf/4Echx2PvdJ1WEaewF2AP35Ki1uA1+64atDaglsFc57O0ILWIWHaqGG25/TwtxfPYDnz+4dPaQgjkCFEIsq34nNUcULDx8nO07HQxWpj2gpjrKkxL8b6rk9R+tFxr8drNAbUPrUQM9n7biVzY7gQEz5HjKpE/g3cgzkvg9HX93p6AOo3XZpwqKKNEJKKnjqG8OUHmpBqvvYMysG9Qp/aiIueh3JYoHr0cm4wWb4j6Gy3HcezAgd6Kj3sEIKWTutJdd3x3U5WCeDD4xFb+YWs0NXxQ6qORikcYnFC6+6uEQRoM1i7Khde48IF5h6kdk7CnIzDOvqNi96AisZdcPKVSMtDEYW+8/Y2Kxv9Z8diZL9SH+D8Tf2GHuo0C+hCV/uVNFYcA86Ml8Z5UXIlrwZmjwrCHz0MQ7fu+sAjNLkMbPqpIky0VxMy0ajQ3rQLNtTwrVuPzI7MQnDfVua8M2oOKWMUtg7WQWJgZRDPTq+e/x+pbFSkkuywd/PU5MR+BvxevdeSorQiKd1C5CUeVAwHDo6BPr6kfbhv09S18P6eKQsB8iNed/1rZ0XLE4rtf+Vj9gv87hi4pM1WdgISDGgzMPR91gDpK5rIHB+oosR/KtqM6O1zrLaGf+gI78ZOZKM+xyo6d/zPU44MByjbtaVA8DYr/CoqMucH0A/5gfC7jwTQ9KMX4OMQZfNB3ixTnENxluPvmwofyYeuI3iqA1/D5C89GZUdgXEFznBhQP+AL+hAEriz13j9ItB6ggq+q7fHoHvD/Fq+KjUEwndQDRsyJQgzHc0K9aisqcDcCZ65FrUZXEPq4JG3+rKKYZKlMC66+QhvvtxCbqQ3NvdyQFviLidWosB0Iqh1xah/6833ohtbivRd4XxAzMl+hzTltaNODFuLtG3VpirMqrX23ogjNKU1r8R3n2dTOc0hqsZ3veFSYCYGM757kVuePGXLyjYeyzIRK3obnOC+oy/cxlEQS7T+nzC49NnP8rLqA7+hBjuHZqgI8bEUdRQQrTy9DDli5eucagXpDN3NviTo/WY9jfD5Z/TMDddMDAdtXXVzAQbEwmEztDcEnt66g0pNGqger9DgEjzAYJ5Vc8dw6MDv6Jg8+OH6wisJrHUsRLNtVCvY4O/bGn7NVIvrXLOwpt56oOOxJ+1BQ3k/VcSSRB1CA0U4VQIB4UHFWt4pA4TkydW9m6/DYbOzN4U+KihAqTYvPV6KZR2qLzY9b0vIfWtOEtGZiiKkJjXE3pknWqmKerTCtR7Dt8qvd6LPXYHvnWZ+lKdaKYoK/sbxuXHIzsfjb1rQ5uxXNO/eyyqiEzy1L4TnFxJyEQrJde73KivwkHC2Fxa08+LxM5BZRDrUPzzmBUYnvc8WtX0SqHkhj59CCmlEWP2UrTsZTEHyOncqMio5C5dqHq9T+2op7BqgUtAwW1FvWTPXjatTzGj1PIudKRulJNHMndXYxgmLDhGC6jISEcVfBwfB4reLD1l50Id9uVzzCSF+IjZ+oODLj0TRdnKkSIA6Ra96lnNgjvCjsO98rK5r+dxEMO+JV0pcoEO4yPlAbUCEHESgH3cqJpjQahRjuVntQuPtQkNvcaiUCY8P1ErT6dmUx3VWR5p1+kTbltBbj05vR6FS0Dr8Gi7dpCE243JWmBOqJuUmlaHYUuotEtRqfNfswWglfCZrqqytGZbanRTSU5t/vI8Yde41Ge5vSxuw2YvHNV2hKfCVaeauy2PxBaVrm1NaCHfTqdoc6VBQ+w4ZKD+NWBA4hyNdhFBJ9Sx1H+aTdQ6sZp26g1UhE93H+fWXH9zZhlGHep+Iw/HTNQWW+qWTn7IPXjFNe7MxP8J3u7lIpqNOMNwsOi2fhOR+jpfl0tNo84WlQPA2K/wqK1CXBMjmGxffXpuXUbHVzjR66zjuXYeOAWKwJJeOuRLqMWSoFeUfMMgTHIRURgiYMfXSUTd3BF45FRQe+Uim/oWIRHOHvKfufqGQkmodS1E50P/vRpO8JqPX4eykKZHummu98niaaq9CqrxuITfcbo/9vSRMym4k3P+lG858MoRFne4muUR2oj7UDDbE0F9MSX6Y3E4JoXqCUmJ9Ukt5yBtE0b20abmsiBno6yeuGHu8p5j4cQrO/70UjHU3ERkJ38qBx/jZMtFWn+Y4XaCt2ALYYO9JafMfdAbUP32c3gv9QquL8xY68KvJjtQKBkHgXw9ef1UEEzgfIxQ5b1enLCASUgeOdXKtRJ5xszlaRExAI6E5s41Uyupf4eQZd2qT43NQnW9DNzFB82j1uQMG5knMrEBSZy5FoIpdgPEHmD0SPZ4biN0h4G3v3IfUwUY9D8AhDIH9wLtXWgcVCBILGtV/ZkAPcuoM9IFG9+xiVi73nWLbiLN+JAon6Qa1BRW/Mwh6boLiAVqAg1iapVbDQ9QzNtz0n5rhL0ZyLTWjx163ERiSSq/56nYZ7W4r52UOol6UjjflinJiYPYsGXxlNw9M6iSHmpjQxuSlN8jdQPkhtSsMtqPDUjqLvkcE0KWcOjbs7UXSP7UBz/xlCY7NeE0vvtpXPXXavtZh/oxm95S5LCxGwwvUsrURruSZJrbDp91oWrzZhJ1mN4An/Rjkf6ZHWo48UB85V7CwJWerOZ4ac8XWHquhlBeXPOKfjk5BpS5SfW32U/ScRKn5B7rmSueoeguoSeog7E9XF1RwUK4Pzp8/9ucKQyTHf71A3tuqJrV+9yrdTD0qde0d5EHnunQXdBbcOHAyhW5QNlc4nsi7fVY5z+H8EifUTdeAKWpaf0GTHq5AHz9GSI8/RUuwdbK6lMM1zFaXZrmJilj+Ipnur0oS4mmLuqYY0yNSURqe1FePSmtMoTxNa8EcfMepyMI28M4pG/vKW6OgcSh0sA6lvSj/Rz92NBga60CB3O+VqRwOTu1A/V1fq7e0lunqGUAe8bvjdN8W4exNp8NHutOjJADE4rhFNOfkqjfC1FoOtzWnu6UbYvlpimrcazfSWp9me4mIufx/zc7QMFcWWZGDo/M/zMrRlYV9jB0HSHX5FZfyIwODW4m9lRisTQHCZ9inPAZQryj1+t3Kh3hI3oxvakwvP8aKOfnCqL6LR2mDH/+OAsqFleYy/4/uoi2ufBsXToPivoMhYFSwTbBnP/eMpdE8iFM+H4GPrPOWM3UelXnehmdqjOAh42OngoaJJcwjuMmwWxSezeBj56WO1HQnXefSPmxOVB4+tQZMYimBga96vSDNc5WhmYlkxI7EmjYmtQ+Osr6hAE5qQ0YpG+tuIgXFtaNKNXsgbholxx1+niVltaUHOINHT3AlN/yzqmDBEDH68gXp8+DZ1Th4ruiWPgWHUzdZX2SEwFEZT1/TxovuNeTToyUZqZxssJubMpm4xnfD+Q8SMW+1o3JG2NB+/syl3+tBAU1sa4VXjj7TGdjelcba6YkzcywiUWgiUIDHDXoZW3qhE4Y+LiCXY+bwoo9VmdYbQ5SbogT05uOdGvSBIHFZlx850GOVu52MfkIruyRlSkIieNeu5EscWlY2dlE+359VxwhQExXrs3EPUxXUIChfGpUeHGYIn2WZyArJAcd/E8yE8uW5iY/wxGDUkqCsXDJncwiMMxkkl5xB58x84IPaHYaMTlPOkQbuQgR8+qTjJ2o7MejH2FjYnuSwt/aAJTc1oLEb6X6WhgTdoyq0+oq+1HfXz9aA+ab1Fr8Te9Nr212nynWAxIvFVGosKGJnZTgwKdKa+xwZRF+8Q0fPjxdQmbTa1SZkrWnnepFYWVPy1xaLHu8updfx4au2eQa2T54iWnhkU/NVKam/qL/qfGk79HJ1o7OnOYnRSMxrmaEpzfxssWq1rRT29vSk4VfVPwvba2tOkq8FiWHoHtCgtaWKgkVj9dQuaag2iFaeKCD4Btwl5xIF0tR9/244jX/OrfajwSASCh3eqRC1vH3ISt0/xJONzp7CzupUHSenneJ+Ejer0XgQV8sSsJco3Vc+VJA9VR5BjGta3gsnc3hABZKEy9WuK4laDJ8f4tyueVp/mLjhs7UvWg1J5Q04eYfAQMq+7+BSjjXBs5Mdfqhj8fzz+f79VhWcgIBCtK3xqTlJZmprWmN6600XwyKF/IprzjL5iyqOZNCEHe+5HU0VHxwDqe2YMmvtOgg8ybctuTsu/UP2imlNvc0fqmT5UtPNPpTZXN1PnbLPKiaVuT6Kp7YWVovXxxdTtURR1weOdc/j/zfT6J/vojeQZ1MU1WPSMxec429Lb7zcX23Ka0UhnYxqc1FH0PTGKOiQMoAHvThLjsb1TnsyigScHiD6OzjqC+bK7mJjcmGb5KtBSbyHFe35WwcQfG3a0CLS4n36rQu1arjd+VnEImgg87/rn6me0MlmomyNp6jSSdytSgIzdKoVHI0go06equP54DPUe10G5Zz4NiqdB8V9BwZfyxTczBF+c4+qByu+jeE6lC90IT6hlXoy5z5zA89IVn8vgQ9e2dMXHIfjgFCeUjHMIPsUdzcECp89jQ48VPH8XAmNZFILB9KyYZK9Go5MxlHw8RPRChY79HgX7eI7oljycOnpGUddz08XQbPT1nrForvuLUZkYajqa0GjkG2zsqS7Uy9WfXo8bKZpZZlB989tU37dGNEpaRw0di6nD/SjR8e8YauxZSo18a6m+f42oGzmfmtlnUTv3eNHLP4hGHulMo5BUshEJ6K5SO1APW7BoZx9Bw7I3U/crM0XHxNHU2TmYxj2cI8b/OhXD5M4058EQMdzfisbH16JZ8c+JVSinHeDjnQ54vkUGgsTqVTws3YFyzZtYfBN52R3smK4TyuLRcyUnzys5V4J6S92lnOg6/KjXd3qqI9wIdERj0ELFj0dQrBwaTI9eM0RU89wTXH0VT7J1LTdkthQLIEfg+RB55zL45FbWexoIjGcc7UUi6ryjOKncgRzC5FXHkT9cuoH3TVO78dgqtCRz7C+IseaXaWjqG9I6sDHfTqRRv82mN2IGin7fLqEBKPDe99aLtqlvUaez86iza4jojr24Z+Qb1MveQ/Q7Npg6JQyi1wPTRYuk5dT4Yji9+meqaJWdRq1+81IDBAGr71pBrX5yUGs8/upfqaLpHTM1w+vaeKaJ9uZB1P/4UOqd0FP0iHiDuoS/Qd28g0Wn08hXkmZQz6/XiAFPttCAn1ZQO8sQMfzuWzQOgdHH1UWMPt2FxlleoVnWYmIVym8XcoWYFHUMO9ela6gbp9qO8rrAO+Pn6gDKMCxJg4PxuRKTQ1sM5kLL7kdLY92ubKjPyEnYSQepxF4IhpYGPXhdedBIGElTgim9qSGy0FIEglGhI1U4goIn2OZNn+PZUjzx4/a3ypmlJ7byDlsno4s4nl0w5OQRRjiifa9F8ckhDoYDNrXdo5Nh5niLiQnOujT+fE/qn9VPjH+IgIgdSIP+WSd6fruWmsdMQlO/QnS+H0ItEufSG55JoodvKHWKH0xtI4aIdokTqFPWLGpmnS3qWZfTSwnbqM7RSFHr8DpqcC1GgwNa56RR4/fiqXb0BnopLVTUtmyjuvbV1Nz7tuh4ZCa6IWzD4aGiQxy6lJRh1M41TjRPeIu6PAij166sEc2iJlL3z1bRIIx8WDvrUBrzJ7qTM4PE+Eu9aaK/Ac1JLCk40eRD+DtjlSVDu4hQl4pAgr4R5WZ6Xx1D+Wf8UHBY/Po/qAeMNs6+q9wImKhQ1N8uZVmFepyFeh6v/Nj5z3Ux8ifZHJ3yNCieBsV/BQUvB5CMPoXFdjLkUr7ABMXH1O14E+t+FXNYp9Dlncvg+RB8+jvpC8XnMnhszQelGB+HWI0vEJauOIfgLmOHRy1DkK1J1QST8SHngUgYJz14U/TIGkV9vlxEPb9ZI171zKauTyKo/a+q6ZE11PL8JnrVPEm0T59BbZwzqFl8Lsdcau5bRI3sC0X9pN308lkXNfj9imicc5PqvZdILzt3iDq+PVT3ig2P36AG96+Kuu+lUV28roF5gWiRvIyaOedRs7gZoo1jBhLRadTaNlG0OrOBmqavpnY/houu2VHUMnk+df14pej/00rq4hlKk7NniQGBrjQ2uTnNTq4sFiUUksP/O9xqD8prJyflWWoluhAvdrxD7yp+ru+uHuRioYl66UO0Q51G18MHEfPmeDo3aL1mTFHHhiG4umj9swQMPAxeHyIdfQvji33dYxEQk1XaIkOuQ7AdVA5OFi8aMluKXftL50PwmU7GJ7f4XEbMz4oPTG3j4xARammU5hBrkxUHxEL3s/SWp6oY6WlOg470ovGP5ogO7hE0MGcLNTdNEl0eHqZXL2yl+u7VovWvLqrnQ25xaZto5l9CTZJWUd3YJSoOiaJ9A9XP3CtqO/ZSNXsEVfNYRXVnFCr+CtX/KEnUvW6nhg+uUg1/HFVNyGUKo9qug3j9PvGKbQPVNa+lutFLRBPfKmruX0ptL28T9VyrqPV9LzVMWi+an9xIXZ5EYfuniP6PNlPHwBgac3+26Jvci4ZYWtBke3WxyIdkM6WgjFahYpdGo5wOqe1oeVeh7CPuKhd2vOVxBedKQlCm53/WBJ9dvKEHufJmzHNQpCxF/Y1XfKW7BXWflsvKk2x4sZC8U+E+RI5pBgJkmeKrtlybC96QZ1zzBNuzXyiZMfVJQVPGZzr5xFbeYWs+Srk4/bn8IedsU2EklEXzu4s5SeUkGOacrC8Gx79Kwb7eNOSz6aL7lbfQQqym1y+tEu1+CqcGiWuQCKaLBhcjqMGZUKrr2iDqJW2nV/y76UWbqhG6lWp70U04Q0V16yGqdvwU1f72rnjp3pdU3REhLYag61QzzUIvffUx1f7+nqhx5QZVM+2nmq5QUdsbRTXCtqNb2SVeSdpL9fzbqV7CRtHg1CGqdyqUWj1OF40yt9Jrnx+gdh9sEl3eW0q9PpxPfS+PFz38vWiguTXNPd1QzHBWo9mB8iif4mJeQlGaGVNYTqqxhUmFKezR87QAXQkLQT3Efo2f5xVfzBSNEeGN35QT3YcXrU1SnOL6dGBYmneppmUmdvxpeN1kFcWTbJ4GxdOg+J+gCFsZTI45hghwswJ80Y7AEMaNboMv0mFyfUJSwXF4Ti55Ch1PjmE8H2J1hp7YYmvfr0hzAuVobqCsmITmkYednFCyiSlNaRS6DA4GNuliV+qZ2Iu6JA4W/e4up9dS3qIO3+0RLc5twZDRjmA4LOpecVGdEyaqbd4npPKTTVTNFi0qm6KoUkwUVY43K6ebgkLNFLQvStT++nuqHnBS/XtnRcMH16hKXBi9+NPPVG5fpAg6aKLKLm/+e1SKiaHKcVFU1apqJcdTLUcY1bLsEy9nRFKDT1KpbtZB0fq+h5pmrqdOf4SIltZJNOghEk5zf9EjsQ9NPteDBplfFSOczWlCSjOa4K4nxsbXoYmWGlKOUpYpZeVcSeijImIFuhg+3b7KpPjalS3obj4ixYlmCp7jiFQ+5IXeHajr7cq6Al0IgoQv6WRxaxAUpq3BsnoM4wVD+KznCZO6zOcz4vSqLZaCZPEEMltziuJJtu/g56EsxVdfbUTiuZTPdMIMZ1la9kETmpam+MDU0JQ3aBxGGGyguxMNOtpLWgfWAwHBcxb6nR0t2ntGU8es2dTCo1qd3UwNfOuo4fnDos7ZeHrxlJOqOs2i+qkzVCnBTRVtqaLcIQ+VN5+i8s4zolx4gCokvkvlw52ixrVbVC0hil655lOXPVQtPoRqXL9BQQeiRcWUj6hcmJ/KO84oG97nkJcqmJMFf171cxeoisMiah1NoFcuJFCDk2GifmATtTi2jl5NXig6YDTUITCBgrNGigm/T6Ue3oIyGIicaqCHz6v0FMPS29OoQCuaFGgkVn3egqbZytOqc0UE52XbUPa7fCr8iEGRqDd7pkpHDnIUzp1QPiSdWdiBr3hU4j6DToWgbjcpx04EBS/G9RhDHha/QRcL+f6s4iu/+SLYvMv4eGkguVDHqfjoWYxLT+0yTnJ4cswSt5rpLSMntvIPWz/Ww9b9j/UTk/9+EwnlXBr6+XQhrcM5BML+9qKXoy91S8RwzzNBtEqYTG0vbkYzvUG86D5EVW1xVPPDL0RQKCrmwk9U+eJDUTYkncrsSqIye46LslHvUan1dgo67BVVE730EkYjlVdNV0snU52LXqriTsD/u0Tp9RYqF/MBld57UpTZm0pl9/qp8rm/RNUrf1BQuJ2qX70tqjhtVMt+CN3GXvFqxkpq5Z6O7Z8kOnuGUU98r3Z7Xhf9zo+kTq5BNOjjaWLcP3No8qOZNPDkQNHH1Znm/jM0/7D4hEATmplUkRZilMJWo/UOx454wKliUe7xqBOnW8UiCJKx8757ReV8h24FwXHWq3IwOklAQ/B3srpkfhoUT4PiP4MiJpguoJ9hvMbUEzwp/qDKWx+CL/RlPLnDY9eLdJgJXQtfzufMtQMbsALdzhxLYTE9saac/p73ZKjgcxljvp0oB6VYz6xR1ME1grpffUtwDtHOPZp6OvuIkce7U4+4TtQzeZholzYNw79FVNe5SdSyHaIa589T0KF4UfX6P2juv6DSG22inPlLCvLnUNm4H0Sp3deo9DtXECh2USHSSlXjI6iGPUzZwqiq5TBVjHdQuXfcovT+c1Rq1xUqG/OdKOfNoXIJd/H+CaKC+zZVvZGNgLSJmjdvU007kl/vVtHUt4TaZcykjrH9RK+EbjT6dA/qaest2jlGUP8fV1L3azNFB+8o6pY6gsY+mC3G/TyFelm70Jx/hohhXnQlCS/SbOvzgi+G3ocd0ns6F+/cyPWun1aZ+DsywqA/f1LZP2lgeGLUfQwWfjxWMD3iog1B4UNL4cAfjBfp4jWmePUYlrc+hNWq+GJfvrbz4wvqEvoqvi6Dz4wynmTLcyrnOouK0TEvy3yIfnymE/jklpzLMA0Sfb5aTAOzt8gIg3FS2fHYHOrqHSZ6xHWmHjGdqF1UX9HGNZ0au5dS3eQ9ogYSvEo2F1W59LOomPw5ld7mpyBftihz+AsqPi+eiq84Lkps/ZSKzXZQqU0eERSdRBVCIyFcVJSfhykoxk+lNjpFsTlWKrHtcyq+8rQoNs9MZQ5+IMHGyuzOpArOd6nq1T9FxTgb1fBE0yspB0RD+1Ikl9Po9cgBItjVk7pHtKfuiUNExyNvURvfDOr+xSoxAOXR7x52DvMQMfz7t2j879Oot62TGHu+K42z1KG57qKCD2btQ8VH+VUC9v5LJw36+bZy2VBXGKHEmdQn7+v6GL98pOyhaAgQGJmhyhmKoDi4JZgeoKtg7yGR5AXH/r2cUFJqwVIAfNU0N0V5l/GlIRD4Ip38KfgBvdh2tkyhK0rjrXVkcky/dDX+yRx6I3ogDf5nnejxzRpqYZqcP+TsyCMM9+z87qJnynAEQz/qemyGaBY9jerFLaWXvPtF9SQLlY9B83/psSi9EV2D7wmV3HxGFF91gkoefEjFtv8uXlh6i4quvEYlFttE6fUmKnsQrcGWw2pzOJXd78RPC5Vc6RLFll+kFxbfpGJbfxYlDz6mEusuUok1qSIoMYdKb3JS5bN/Ck5gq7nMGA3tEXVNK6iZaSZ1OzNbvBHZl3qljUALMVa0SHiTOt7dR69fWyuaozy6/Z/D4kPQYsyhwZeGiIlXetM4V32alVhKrEX9bEcg7ItVcWjJzfj540eKLz7ievvxoeKTkh50Lw/vqg8v4jHs6E+uKPvBp0HxNCj+Kyh4iWBewo+Z0Xxk3/u/a0z9iHHv+Q+UHd2DD4ESF6c42YxA/7UnRvF1GTwNf5YvSPCcSp5Cx5NjWLfAMDn93fObtYKHmXzomg9KsRbnt1ArPhaRMEW0S5tOrdFlNItRbS9spLr2tVQzaruommCmIFRquYMpyvYVlYn8loovPypKHMimIjOzqPBYl1rwFT07yk1F5/lFyY1HqOTaJCq5JrHAWnQdW4+jm3CKwqMdeN3XVHi8TxSZlkIl9j2i4msvitL7b1OQ51d0Ix5RLjyFKsfGUs2IHaIecp9WWavzv0Nr53R6IwU/7RNFqxPrqcWZLfTGdyGi6+NIapm8gLp9vFIM+BU7S8Lg/HMlfZ2daVzyq/Smt5JY4ilEKxAEe8wqnPM8BErEYcVT9XgZyKzLipda+gzBkLdG1xPUtwt19+C6upTIQeEMpnfxJuwWslJep5JnVLEfEBC8xhSvHiMryPAJMLQGfOU3O4yA4Gs7+aottgwJzduuQjTdW0VMyGxF/fzdZbYU4wkyA7M3U/PoiYJPbvG5DD5KyfjAVP3EddTm4hbR3L8YOcQyqm9aKura11E9/zaqZdovKoWHUnn7GSSNfsF7bbFZUVQyJFsUeesIFZ73CRVa/EQY466RMfk6PT/OJIrOSaUXZqZRkamOAjPxGB4vMj5WFJp0kYwxF6nQokfi2bfv0vOTvVTiYI4oNjOagpB8lt7uFxUzPqPK8XFU2xMi6jgxUmJxy0Uj1zJqhu/V9upWUd+F0ckvTmqYslHwuZKuOdHUNHKi4HMlHfxjaNzjeaJvai8andKGJtlqiDlxhWl5NEZ+ScpxxKBQBIY/Q8Uk6sIppz9ThxE43z3C75fUx2jxf/oQPUC0uhRAULjigunQZkPwCrHuhIIlCc9hGJrGwfBQmZHAuNDN5C0FEJaok2q4+RIIEL5IZ4KppuAZ1zzBdvCHUwXPlup1dx21Pb9C8JnO+u5VMqmF8VFKPihVz7FBYS97JXkvuooDQlqI+P1ULWaHqOK0U8XA+2juk0Q5x99UfClaid3/iMKjrPQMB0Nvixr2ARm9AlRoZJIa6qRCw9LJGJKcj/8uNNRFhUYFhNHLQ8bIz8noHiMKLXxIhSemoSv5RRRfcwGjku+pzDsXRAXXOaoUb6VqUTtFLdtBqhGxnV5y7hOvpOyjuu7NelINGp4KpfqnQwoOi6dvode/Csk/LN71znLq8/nb1PP4KNHD1Utmi4+zvCT4yvjtCIStbrWXDw+gLhxZKhEDhmi07pf+ULzcQQADhJ+fKBOST16xL2/EmZ7wNCieBsV/BYXPGkwJaPYZ5xO8tDEvDsrivXpu3nRM8RpTvKTQAZ/iySDrERSr3Iqv7ZzrLUVzTzYUA2PbyDT8jo7+YuiTjdQ25U3qcj9E8HyI1r86qf6Fw6LuVZceunaHiFq2EBl2Vk+yiqoJFnQZYVTZYRdBu/ZT0OEABSV8K0rtOEsltn9BLyy8IJ6b/wUZE99FZeca9BMZXY+S0cWuOiFQ2iVA4r84yOho0f9j3Y6Q0f87MkZ9q0aepsJLf6ciM1JFyQN/UMkVfiofuC/KveOloH2hCNgEUSkMv9vMVN1vFTWch6mWHV2LQ718IoYafJxGdY8dFG3+8VMj33Lq9iRSvGqeSENpC7W3DRTdfUjazW1pzulG4k17aVrgKkLLUQ+Ml0HYglzvHb8KO4pE8jf8TFbYLaR74dSA8YL1j5Bj8PU7LMuPoEh1BMtK+OwMIorXus44plx+PX6ecEbxomMh+LnOp/bg7znWQnLlN5tkqUJzLzWhwfFNxaSbveS6jH5nRov2nrEyp7KFb65odX4T1fevpwanQ0WdkyZ68ZSLqtjjBB+Y4uMQFWJUOSSVnENw68CCwpMwekjBqMMhyrxzmYrNtlHR5ZdE4VF25ASXyAhOVl2zyGiRREbr26ppBhkNEAT1fLn8+NuJx/HcNh+oVgii9l68/ogoNP4CFR6TSEUXnxPF5mL0suOkjEBY2XfQYoUmUQXfLeW6KNtdIVpVMjuoxoXLkiSzWhlWeuWam+pnHRQN07dT84y1SDYXiXZpb1F7/wQKPjlKdIrqSiOPdZE5GGz2xaY00Vqd1t0sLea7nqENyCU2Z6pdyP3saPFDXMp/AS2HxaAT55UTueB1tP6WQ+p8KoIiYnsw2bcawhdpyMLnfo/iF/BKtbwED+MV6LYi2lw5ai6y1/AnL9Dq25XEqq8ayIW+eZfx8VVbk+8EU39XJ9HL05+6uIbkH/JtaZlEbS9upXp82htqW/ZSVYe54LA1H6W89AtVwXCTlTuUggwfCd3mgCjHrcMGJ1W5+JsoewAjh00uKrHMKYpvvEWFh0ajMq2qPVqCJqjkWl714vtkvPw3GbWfFHj5Hzz+CRk1E1Xj42S8hkDpGS+eHW6l4pvvUPEFZlFqs5fK7nVTlfM/ilLrLRTkvkelt6aKMjv9elj8wj+iCh8WD7HkHxav6nXI2dX6iTtEq9MbqHncROxAk0Un1zDq7RlA/VydxYQbPWnuY5TrtY5i9s0WtOGPJvRmYkURmlOSlqCObI8VHwbfiZTAkqx4Fr7Nhvr1qgQ0Bn6w71SxfELsaVA8DYr/CQq+xVLcIkOk70EesV2fxM6iC4lDP8OLljJeL5KXJeQFxxivMcVLCk13VhSb/mgsl+WPT22ujr9OwxNb0ND4JmJUZkfqEdtJJtiy9unTMexcIsNMGWo6Q6j6qdMoNLOQcxl86DrvXIbtKwry5SCh/EeU2nFOuoyy+32iUvI1KrPDQqX3nBBFF2YgWUS30M2lmtjIqIGKrp6jKnxGRgkzGS/E58LvJUxkBN0ho1p2LiSbTZB3dPOLQqMDVHReGpXamCRKrz1MldJuUlBkpii97xyV2nKSytn/EjxcDXLcze9eyrvf13MlYXZR8/ptquVBPuXeKpoFlstxjE7WAaKHpRuNSu9Eg+KaiuGO5jTlagcabm8iNuS8TltyWtDkhGri4OOKtPn2M7T9oiE2cjcSi53XoqwIkDPIM/y88wNflMzXg8SvULZtCAq+51baQkOETsVPBEZqjOLFyHkNSF7FlvFCovuvGjTPrkKflKTF5yrJ6jGMFwzh9SFGeXjyDHKKY23l2k6+aosNszeW6zI6xQ8SbZzTqWnSSnolabeoFTBRpQTkFBd+FnJyS85lPBFlIr9BzhClIwzgpLLYLBtaC7eQgNjmw3NixXOj0Ur0T9CcgDVCQFQ6QUbhWGVcA7QMRs6/PATkG8/GqwqnkGvgda85VC8TFR4eTUU5d4EyezPxuSYESIIoNttCJXZ+Q8WXHRPFZkZQmdCP5RiKHP3kcyWOW1Tl6p+iotVLNf3x9Ip/n2hgW0Yt7TOoXdQg0S2qA41wvEpbUX5stKcRzX6/o1wNx1b81k4WTlnxaX0xJ6WSLM42D3XFwj40aAd+RvqUDb8nYdQYiFPWjQiWtww6vUpF8LWkfAO2vNsrmeejZVhv5N8agVfC59sh5C05zKMNXngjb0nCzZ9XoJmZtfOXE5qQ0YxGJLbKXwqAr/weeaQ9Lfu8hRjjb029bT3yp+A3N82gujGL86fPVbVFUUVbSv6p79IbbHJiq+Tm06L48iNyUCpvyMkjjKLLLyNAHKL0npNUbE4sFX07SzwzFnt5V7QOTU2qGkYTz/v+FQA/AFoAw/IvGIIaXxc851kkppXRgjRNUF3deN90emFWQBSbjYDYlkYl1yWJokvPUpF52NadD0TJQzlUfDVGRatTBLd0fHa18uk/RPloP1UxRaIL2S3qxi2nZqZZ1DZ2pAgO9KORyW1p/c+txGBHaynXhX/3E6NS2/yf5ZYWXKhDq69VpJD7z4sVGCzsQN1Fpii+wwDPoIs6oEybUMdvo75nqKPbnwbF06D4r6CI24juY5khApvRnCDh5JunMF4UlG9jwAufM17rmpc2DssuLSaaK9Omxy1pQnpzIWtMZQ+h0Vd6C14fYnCgM/WPaSHGneoql/LxRTqsmWMevWJaLRNsWZW4KJlCVy4kXZm/0NPfK4+LEu/wuYwjVHi0Tc3/kgqPtEtCyTiH4C7jmTEBYfRBIHThnCBVVUQXYpwEey4MM40HQP/C3clpiM91g4zy3PUkqw74vWcUPTM+XRQeGUVF56fh82+IwiMtVHjhdxi2ukSRaQE5V1Jiw1VRet91CnL9ggTULYIOp1Kl2GiqZd4v6rk3UzP3AmoXmCL6Hh1Ko452ob6RzcTwzC7U09SJRpztKeY+yl2Dy9lU8M45xY3c4lE5sRA79QZ0F7xGKYtCThGBXCLqHeXca1DWLp2fy5xbERSHkWgm7zCECQHBd9bxOlRYmCaYW3JxJstrXW/4qKJYcOUVWv5DGxqVouSiYEtHGnlntOAFQ3h9iN54TLgHyLWdnbJmiuaJi+S6jBe9EaJSdBSVt5ySjJ3xXlV8rolKHHwoiszKosJzP5GjlHKkcuK7VGjiJRlhsEKj0ff3R7/f1a66YC9v6cZIw6KKcSuBv40/c3FAYHRhOHOhJZBAuQ+PcuH5RdB6VItRr+L3TkkFB8D62pHMJtGzw2JEoQlnyRh3hQoteiwK87mSKT5sf7YoPj+BynmQH+1JFxV816lynIlqxu4XdV0ICv8iubCJdbAOpt6+gdTH1ln08vVBuc6mMV+OFz1MHeVC5Qln24mFn7eiJR/UpyVnK4j5gedoUTQSzggVbkKLgUHEkSTFAwrTNq1/iYEdCIowDEHyMlFLiCaXPHOH8XCUb56yOUItRJAsSCmSvyTxxvsYYaQ1o7m/BAuO3NFfjMtfTqiTYyh18SBQ0oeJN+JHyWX9zW1zRCP7ArlIh0cdjGdLV3CeQfJ2XPBsqeLL0Urs+F0UHuukQouzCw5bD34XP7H39rEqHmFwQtk0PlcKhplo+qtkqkK856PlMD7OhZZEupC87oRbil9JW5GbudIA71H5iKoZn9vyWNUbCKZuCMbe8aoXgmj4R2T0iFUL/6HnJh+hYhu+FyXWXKYyYZ9j1HReVHBfoUp2h4xAWL2U3VQ/Zi419ywUryXOoDcSxlF3/xDRLTCEOtoLllsa//NkGpTenRagK2GjXE1py5OWNMVeWfBi84uidLEYxks2ph0tONOdgr/jQgoGF5G7nwbF06D4r6CIiwgmB8avzIv+JwFBkLdudCT6H77FEt88hS1F17LsVFFZ+JxtwVBoCMbOb9NQwcOnSdmzqINrqBj0aAP1+ngRtUuaLppbZ1CLwAqqZ1V8bSdfylfDelBUdnrkXEbZ6Nui9O5rVGLrJ/TCkluC50MY46/piS0m5zKyEAg+1RRJZSP8rGZWFRO1y+BgkIBA829gSGo4cnEg3ILUXAgi43Lu48m5rsAhvN6sXvDp+1a1qgb8uXivN9yKz630/5aM0ffUyNP07MIf6IU5p0SJ3T9QyfUZcuyC8bmSyp4kqm4JFS96Q6iefxfVs6wSzX0rsAPNRH4xQ/R8f5GuwZU7XW9SzlzqGtuJFuYMFUPszZBwtqWZGTXFmmslaS0GCAdSFK8/HoMd3+pRvH4WL+2cd2M8SxSCYveu4Pz1klKzEDX8+2nFN2ELO2nQphTFC5xvvf0czU2voC41oDGexjTxchchC5deGU0dzQNFzw/fljWm2l7dIhpY3pb1IV5K2C7qnHPLtZ3Vj58W5UPMcl0Gz7hmpQ9exmgCY/+V18Szo91UaPJ17JFJiiughV+PUrIaXGGo9CI+ZZwizSG4dWCcZOYdn2AZcAx+y8WtxCXghPRurj1wJvf1jEcnyBueS1IV0HrUwDY0ylSt8Pkd8bxeaWrcGXp2rB9J8HFRbIGfSu08RaW3JonyCVeo/GEnVc3MUhiB1bngoRcTdorG16KovmUJtb64VbR06Rpcb8T0F0PfHU/B1k406Xp3McLWmOZfbURvuiuKXR8XoeU+vcMRO5Bm0CGIDSi+fkfujBBQh8MRFDNnBuff5pnv6ss3ceVVcdnRPw3anWHIHXTkLjof6zJ/8y3PiIX+ojTZWoWmJtUTQ6wYkqZ2yl+SsHPyGGqTMoc6Z8cLXj2GFwvJWwqg/u9XqJrHQrW//V7IFVnhDgqK8IqyexKo9JZEKrHULoot8MnkGE7smNEZzXyb2wWHrfkopRyUyusOuBv4Cz7KxckkB0NeYpk3HLXmyhuWckDljVC4JcFn8ShEIMGVwMpWz+A5VfFYzYBq+zmCFd3YUPXcWAsVm++j4ottohS+T5mddgqKTFSHLFRufzTV/PRLUeNoMtX7/iy9cskqmv8SoAbpm6j971GiZeocObHYL32QGJ7SmYbaXqXJvvrKWp0W+IrLovNsOVqFPbcxqvhbbUGAnHyM0aRbfUo68eZ3Uhs2Pg0KehoU/xEUe/cE0/lrhrj6iUG+S4bcxJVFowtJ/U2XBWbWXwzacUxvf8DW+nhSTWGa5y8lpiW+JOc4eNFSxutUtva8KQuOMV5jqvWTNFkshPFFvdWdkXKhL6v9zfdyKV9Vf6KoEGWj8jEBKr3BJEptOiJT5WRyDMh8h6bpemKL8bkMqXROHlnecQiuxLxg4S6DK51xohkCO3PtyP0buYmRmIu7C36dK9fPpJ9xJ9cXZJS7QMbLP6iGXpmjUWhYmnhhZqrM+yy9wSKColOofKSFqni9osaVK1Tri2/p5R9VVXMoNcq+SbXiNomWj9KoYfJ6lJ9JtEZQdD06gfrw0pEw2Nyc3kytS3MTS4t51udoNXbcdQ61BTu1/XddJoIl/oD8EUGS+oE6+Zkht/T+8ie1ntfRTE8PJh/6Enb3EVoIj97Vl+1BxfN9O/mOfOwd5BeH8QZ8Nx3GazPxgRFeCZ/x4uey8Eigi+CFS1taxssKdKxRwiJq/buHGl6NFvXe81KD+5epuiNSBZxULSFSFhZhVc18TUYElTvoEXzE8IVZKOzh6cJon4BKgDp/Kz65JQef8io97zjEvVw8mkBiyqMJwQHAgbA81wrYCHup4OAVJ5rnSN+LcULKrUVeS4LfiyLgXnyi+DwJT9YZ5BdFpqCVW+mg0ptjRPmDkcgbYumlM05RzRZONTJ9um4G1P/5AjX4KpPqZR0UbR4lU/24+dQ1O1q0ck6lroFR1NfdVYzxNKOJcRVoTuwzYmOSHnk+nKGizmKwcBx5w4/q2AODdlp0nXRmxvO/QQuSmKyOZyEoMtKCKR6jCvbTQ42avNs8X/0LgWFFoJxTfEe+fTykual2xOpFrXm3RpifXFIWPs9b1ribta8sWpq3JGHHP6KoAbqQvOWEePWY+h8n5S8FUP8HNJvX/VR51TRR0x4mF+iU2XJYlFybKJNrjSEpigu/Pir/xceKz3RKs55X6TzE5G4irzvhCuWEMq+74NaBg2HivywjDZToXJxscouR9x78Ov6Z9xn4/Xk8VitH1UGQvoHWYnBAPD8pgYovclCZjfGiYrSNajgRANsWipfP2Kned6eoEb4/a/L4IlqIDdJCsMZnd9Ib34VThy93qIxZ1C1xMA3wdRKTM5rSLG85WoCdkm3FjrwVQ8v09xRPzzuZg2D4Q203GXQBO/qt+yrSprPseHVkdvH006Cgp0Hx/xMU5jhD/ICm5fgVQ24Ez/h2RbyGUt4tGk8hGd2LIDj7teIxL99iiW+ewvh+GZOQ7OQHhb0v9Xx3ObU6tlh0fBBDDdwr8hcdqePfS3WvJ8hiIazhg+tU54qXKi2dLGraQnODQmlQJBRMtOWg4Cl0eU03z4eQoMhLGDko+NA1/2QcFDzszAsKrnzuMv4dFEtgO0Tl+p7+Oyg4H2H4/bn/Nyjw/MHJ4vnJDgSF83+Couq6OeKlk/GyxFIjuiGaPLpItU0bqeWTdNH49A5q9xOSzI+2ik5Z86ibh4Oio5ic2YRm+wuCYicGCxvDUD/fKr7O9BSC4ujPapdF7+t+83cVg7wDJZJ/wFKCIh1B4UWFsnvoY6ITDHo/Wx3CY6l3DbmrLwu/jP+/qB/E9qdowsl302EznUFyv4yByZ1Ft6Qh1Cp+nKxiy/g6UF6nstF7ZsGLjvEaUzXTrKJKXDhVNR2iOpcSRVXzYbm2k6/aYiXX+SVxkxnXPAu7PVqNhg6dPcVKxJHOh+DkkHF+kUkaCIwPTF2igr38EGyggpxiae7fnFPkPQcJrBzL4PdlfKyCA4+TUZadm1Nkq3qcaHJO4RPPT7ZRiRUYRW2OFeVDY6iqNZrqXPaJ6vZwyZ1q+GIF51i8/laDc4dFm3/wu2WhLKjGWjmmIKcYiUS+qxjrbUaTTBVoVmQhsRF1ths7cN591iIwWAi9gJb9vuIckReku/1E8VXq3/2D+k9Sp04gKEIOBdMpVDR790v8BxKT03dVWKZBGWg9lplUzOfIZtPRKlgUL7g+38pD0xJiqrc2Dbfy6EM3uFvyaIw+ZiBrjhGNZGmiVKoVtV7wKnQ1/LGynBDjFWRqXr9BVdwOUTHeidEHT983C75yS0Yfw1zC6IRKa4JKeukTVZ5HA++Snvpmp0lPbuUlhbyHo9LyE1F+7ABpi5FnH2lrkPccbmk4scz7+zvSQOHhKfuWjDIn0EL8ohq6ZEZ43uijyFvJVHJjKpVeHS2C4pKpQqydqnpcosalS/TivZ/opbtfiGqWUJTLLZTPBtEKo7VGKf8afSRj9HFkPPX19xB8h6Ppvjo0319S2QvTGgTGGpvaiACJ/65gBBnAgMGM3iDjU3UEI5Bj2Nk/+VZt2YygeOutYPoFow4Wh6h5H0GQt8RvJoaj+zDKMH+v+F4UfBebBbZCYqGvKE21Vsy/3xbfXmlkWicKTuwl+NYIbfAlOufEi/pJa2QV27w1KnlJwqoJsfTi9/dEEB+n4BVkcteGKMtLAWxGC7HKJYrNd9Hz4+MKrsnobNEZ17USFU+feza3Sc/v9/lMZ955DD5sfY8KgoSPU/BzYnJxd5F3nIJbHcY5Bbc2Z3NxK8G5BgcKPOMmowo+o0ayaoPg7MJBkSieH2emovMT8+eNlt7iobJ7HRQU6VH7I6nc7giq9eU3onqql+rdO0+vnDWJFr9i9IGg6PBntODjFB0Sx1Gwv78YldGFhlmaofzri8mWKnL8aJ7tGcH19c7HBiXcV9swwjjxSHsB9mmOXumXd5xiDd/v42lQPA2K/wmKHTuC868dTUUg8ELqfBUYC8E4N/Qkstkjiu/tveFyYbmrL+MbuU5Iakyjj7QXvIZVvyODqat7iOAbqLTks3yf7BV1o+ZTsw/M9KJ5m6h7O03WqeQV6FjQQSRiyR/mH+jhSbC8PgRf+c0Kj07Q6zh6exVfqNP6qM64ZtURKBVPkVE4SRnXSQ848ZlOxscc+FwGH6VknHfwYxwojJNKBJl0GRwMbBNpHpL3Hpx03tAZWax8BoIRCWzDTPU6gq8TgqN3QBQad4IKj/FQ0YVHRbG5Liq98xiV3uQSFTy3qHyEm6qfvSCqxh+muu+mUm3rdtHsMzu9ErOQXr+zVzS3TKWeHy7JP7804MwI6pfQicaeUmNcjWkmrzboqyC2597bfcdJtS9N7xIQ7Vd8XU8S6jmQrkJDEBSRocH5Z8h4QRIP3oDvmsvi8Ps+nqDLT4YVZkNOxc7MqC74zr48eZRv0ci6RneQCSAdHEPEoCcbKPjrldQuZYZobp9FzfzLqZ5tlUraRS+6DsmShIxXoeNFx8rGvC94BRleMISXAmB89bcx9qJexsf4bCRP23/NpfhsJQ9Rq5hVBfxfET4h9u/EkvMMDhjGrQkfmMrLF7jCOankVoBbB8YBsZ/yD2Y9jyFueTy3slXVxWua4LG2uTrg/4LRWoz5RQ0/Rs++/ZNMEGIl9vxEJXkJA9dPouwBj5wlrWYKES+6Q6mubxfVs64VTROXUXPrXHojabro/elyGsR3L8hdAH4SYTQS15kW5AwTQxOa04bHbWlmek2x5npJWufVm8mw3cgzTBhRxjkV30gmgeudf4f4SG4p1gfT4U2GiD6E1gKtQWys4msE+NLBzVFqYSikvUBTLBUE3+Z5XHJTmv9HHzHseA9Z3W7492+K9jZseGw/6oohFGvnQnfimUotPG+Lhub5soptLVeokPkUjtNUZu8JFfMdleDVY3LXhuClAPjK77zrOuUyvmBUXE+T6p6kM66b2lWjFLQeGBlUzlR86ltGD3lJIrcifJQyr7thPMLgAOFgYekaDJXSFM++apikAcjaIDC643N6mFRv/D7mZ7RiUWrefXpu2hkquu4rUXztNSod8hGVOXRZVHBcpMoOF9VKPCx4R6kXPZda+JeI1glTqb1nInV2DBYdzP2ovWUoDfl8mpj4x1QakNKNFj4cIPhWE1tyWtFUV2Wx/GwxWszHk5wqDiNGH8TGqVQESmwIGoE9ateGp0HxNCj+KyiikVOk7jBE/E70L9FG/g3nww4XXATEeF3tua5nafXNioJvBr/k29Y07vhrYt7DIdQjriONvzdRcFfS7+Qw6hXbSfT2D8SXGizrPLEWgaX0im29rIrLKsVEyzqVZfakiqDEbPTBZipx6IkoMi1Z5jzyld+MJ7AY4y/QMyOsgi/UMXqbCi7+6YimvSXyhlomVZQrm/OI3CTRQPMuB7R4CJqHg4ZzC04oGd7jOVR+1RjVEoHRHu/RNVfveEl6+coxVmg8tmk48psFf4nCi39FsmnPX8+i+HwHlXM9xvfLEBX9N6iyKZ5qxu0TdT1bqHlgsSwCL6IGUXDKEOpt6iKC0wbQxOzZNOH3KaJ7TEea9cdgmnylg5h3pxUt+7ghLTxaQcxPfo4WRmBoeliFoz5jUJcn0hVfFRjPE3d3KsseBIVlSzC53zbEUTxo32ZQ1D4Vgwjj2dz7E9R2RNhSRFbIk7JiirMq8opWMsWc8axivhXSkKPdRbdoJJ7OjnJ7JTbyaGfqfxyth3+i4NXw65rXUM2wbaJyHEYfh/iUuU8E8Sp0B9/X5YTYvseyYEjhielqyW/07OhEKr75Q8GzqvlCnWfGpQuedS2TbJukKJ4xJSOI2FyciHJw/Lul4EDgZDMmF1qUIARTg4DigOjBn5Ehnh0eQy/MS6USWz8WhUc5qPCyP6nw1BPi+XE2Kr77byqx5bYotf0CBXnvU5mtCSIoMoMqxUZRLcsBwav+NfMspHZJU8SA0yNoeEonGmhpLQb5OlH32M40ILmrmP1gCM37pQ8NNTcRm5DnTfNWp3f+CRILMKJc78FAIV5FIkcMR3BEHVTuAwZl7TXIsVB5eYq/G93HwRGGMC9AU7IegRGu+FpDe+K/rvvw6PJFB395Tqy9XYHmnX2ZNma3Fnyb57FHX6NFT/oLvkUjZ8Rvv99M8O2Vejt6UtvDQwTfGqFu9GKqbdkpeFnjCpbk/DUqZQ0rf07++lIl1l7Q1WO2/SJ4KYCiS85j74sXpTYlYbRipxdmJwmegm904es++BQ7VIXneBia1wrwwS4emuYlojyc5b+5a8ltTZ7xa0KZd91HFw8VGptBL7zpFkVn25E4uqjEYosouuQ0FZmeTMU23RMl3sH2b75FxRfxOZAEvUpsayJVOvGLKB/l+7/XffAaWfGz6XXzaNEnpT8ND7xGK79qKYYFXpdbXS7KHiB4BjcfClj/sI1Y8m5dWny2IoU+eEHwfAq+7iMiF1/3YbcVXMZh2aLXffhn5+L5FE+D4mlQ/E9Q8LWDHvzBDk8z5GLTQLTimd08gTfCrXbhDQ9c0zvtyt12c0rRvMxKtPhmHbGJWsm9vwebGovpV9vR6KSmckc+NiKhidxzq2PcINHaMZ2a+FZS3cAeIdeS2l1U5eofghcuLb0ns+A6zAPvyZpWxddcFCXfuU9F57gRDF5Ret1hKr3vCJ5jEs+NjCajj11nXLOGqOCKR3S2FOMJMhIc/2/38Z0elGJBmTrsbG1VPfla0hgqhmBgZQ6dQvBGU6mNLlFsjo1KHvyLiq+7KorNjqNSOy/mf4eyYeg+TBeo6vWHoqIlkWr4TFTHt080MC+hlrYZ9HrkINE9uiONcr1KW1F+bLS7ES34uisNjGosFnyu91pf/V1DMd1RkUJyytLbyAlZ6B3szKi3qCRldemVf3x/WWbZiHRghkFnV6pIvpY0ZlMwxaCFYOnoW2LxpKMWdQyJiAk/w8xqc6zOttp7U609x1edl6HJ8RXF5r/0qvMpp14VY4+2peGu5mhBGotRaR3kJmy9UoaJdoFpcuMWuTUCyO0Xzp2n8mF2UfVmtixcykc1Ga9Cx1dxl429K3gFGV4wpMwer6iUeovKbOfL+FJF0QWp9MzoJJ2Cz5oiMKr7dfoc4xlTxZBbPG/JhUovGqPnMvgoJavG+Qhe1z1VcFL5wtxUKrnaKUpviKZKmR9QuZBkUXrncSqx1ENlInmlvq+wzTlU3odWb6tXBMWd16vOIzyixrWbGHVEUV33NtEsZTm9jnLJu+Vlt4iOchfFYdYmYrCpEY0/9jqNS20u1tPrtIVa0FR3dXHgYUXa9tEztPGEITZnIZcwYYe3KDtyxJOo17Q4xasM8GoDeSsPJPBV5/HINp1ILhnfqtCJbNTjUAloZo7ijQ/bVYgfH4JE05mj5uH/w3OK0qpblcXqrxvQsnutabhPcTPHd/UdnNRJ9LT3kdO+7TE0Za1sE+m1S1vz19rm9SmqOMz5azcEhVnRYmhXwsqi4nkFujIHL4qgwH29LvP8T0KWAsDffDaVldh4g57l9Sl6mdQbSBIbZxbMp6zzIxkvPSk47c34TGedXwsOW/NRyjYcFDzkNSGxNFNxTioXWwUflSwXEsi/irz0ehuVT/4TLcJ1UXqrj8rsdFKlU78Lbh3KR3mo2tnromqiW9bxqoeA+P/aO+/HqKus/w+7grrroywgYG+LVOkiQoBICwNI7wQSiigkVOm9995CCikzk+ktmfRAQkgCBGmiXxUEG64CyuqCrLRwvu9zTibZfR7/hPzwImQy85n7ued8zj333nPPYTocX0NtEsMpxB8pcOVlzuvBq8VM9M+DaMaDYTSuuLMw80I7Wn37DZpkry/spydpAcfV3lWWwqJvhiwPOhU+y+PA736P4kioyEuyQbHvrFaKaqX4I6VwJxspE/4DE3Dq4lXwDEBalkHqaztLlPgjukHG9SWYzaVQDEcNWvv5k0Kk9Tmaea4VDbG2FrgYPNf+5qq+TGf7CKnb2doxVWhXsoqauJZQs6N7hddyMTU7ZKVnXAonP6+fYpVzEUzdA5kw1V/B1+CTVcfpb7sxfd3Ngb3/sVcSZaG/zC8VHhll0bOdfGqL6XFIz2VwxDXT3EuGJl4NjAnC8RC8/f3WJaXjOV1K50UpUCO8RBKuPj63WGDH9qnNR/D9NuFvOz3Srnq2IqHBoW9k2lk3xinUT+ZEJR/Ts16X8EJ6Cr1WYqYmBXsU7wpqk7eM2niihZDs96mTezQZC0YInXZ2odFFPWm4t70QVdaKIhwv0opP/ibMCfyZlsGh3FamcPR20kn4hGmKNVczCQT3Onz4fxHk7GSDAAq8UIpcj5H8mL8ynB3PZq/Kjsee6hewCImHlLyr+HkCs48MZS2UaBUuutCpzHI+SlNdtWlGcTNhUMpbFObtQ13s7wrDytdIIVcu0ci8kcNF2DzUuGifwBV1/l5mp5dtewXOIsfJz59zmwXOZNvQnEwNzFaBs9DViUmXdELMk+uL6Ikt39NjM4qFmnN+kMO+coyP4b2Srnj6e7CCuKsy43H4XJBghrx3KuDNLeOX/7GXUUiPLLglC2kMO5XsQ7B1EAuxCwqxN0nSJjINLSZqmIh7cJqFF3xJ9JJ9v5R/YP5ebKZmPxRSo5wdQod7AWruWUA9yg8KbZLH0VBaj/4bKPT096OBzo4UXdJcmGiqI7vVcyE3hgOpV0EWu9MVllfRHQPt9SmfQa774Btef6BkZGl2vDy7cjyPlcJtlBrmDOfRtMEZ+eGWcuSsgQ5/ZqAzvysJ2TBJsBZcj5TZA4XYgRnKOpeygbfVnbVofPKLwghvB+qb1YcGnI8Uepz+QKr6djy1VOh8/QA1TVtZmUOSUwY2LY6h5vk7BC6NwJnwg3muOa0xZ7F97uBmgVMS1k87T0/BqWQ4nZBkkNlyW6g5xiUnvw09kxQ+xhfGRwTcyrBsMgz2V8ZTagidvyJAxitw3k1Rhu4JCq9STiqqnHLyDIMdyuBwwdahIefRTNkhcHtfjMeU27RVaJS5m173baCmmVuEZrjfRrm7qN3dbOGNkm301qWdFHJ5vRD6yTzq/8N86pkzTHFyHs2OlXk052Y8RuvzoAgOhcMmOQrfnK2kFsJCQPClPyunr8MAQHb/LFcSEzWPpiteKcusVopqpfgjpeDc3B9h3GEuHjfQta8wVclXrsG8cOH6o5eVz/C7Fa9zcTmG19E5/mKPV5mLMWm2809Vubnz21P/9F407uEMoYt/NA14sJZaxY8XupfHSy7q5vkbBB5KuOYWV9NheLrK9TJ4QYdpwhl4nXBIrXsETlzKeSor81lwbu6psfTE3odCrcq9kjsC75XUCC+FspgFzsPN5zI4uDbIoxOtEkJXa6xZqDGuiAwjCqnGzFvCn+f+SjXDU3VRCvA6RB0PL0j5hQaHvpUh4yUoA/OaeakkImlsXiw0tcynlv751O7wMqGJdzm9edtLb5RuUzIW0zsPD1Kr5InCu3fXUpfAWBp7d4bAublHpnegSMuLQlRiTVoIwcblKRYoQlwqHlr2BUEy5PQ5XIBT1xROrHrlttYrZc5/ZKB/QeaeWEVyc/9nFn9brPoVGWnKhSua9T2YAT4BGshFRnhPhEnDF8bg57YKlkNJFkO5pvrqCbyD2s8aQhMfTBW6OQfRgJ8WSekCpl1gphRy7fjNbqFV/kpqU7BcyitJiaWcD6SASqvkqQJnw+f0zC/GbxCetZklcenftvuFuo4f6an9X9Jfl5QKT+y8R7Um+anmxHyBT3//eUwaPT7TL8gMZbGTnvjQUcUiB/3Pmiw56sc8MsYr8RCPvFci1Bpv072MteeFJzdhFuS/TbU3OoS68TniQ7wA68A0dq+h1ulzqZV5mtDWPJk6ZX0gpa+YtoeWUsvsFfTWl9sFPvDzVvEiCv14vjDo38upi2WQnOZn+rnfofDMtjTF01CY66ohCdu3mRTO4p8D5UhIUJyQxz68FjiisEJ8/6uBTCkKZ/Fnf/LmKaXMzyuaW410B84IUwZH5EyhOh4MT1m4PsS1O8oNaFxOgYEun1HS3dBMEONQuN4H58SKdv9FCLc2osgzGD6ODBDC70ZLYMjg+yuFHhcXU+uUSZV5qENv7qO2gVnU2R0pdEnsT53iB1TWyuDSCI2TF9Crzu0Cp2Z+OsFdmaPyqTUOyX4TzC/1xOrTkj3mL6t/FB6bfpQem3UEwrcKPGOpvS4R01hTFWuT6MlliZXhc4/POSTBMY+vvCLIbuf6T/972XqtF9PNX4W6+7X64UsYKpjXYR1amqZRj7JogasDhVoGUmfPOIHLW4TejqGQi2uFtrYpFHp2vgTSSDCNfRiN+VcUDTs/QhhT3JvCPU1pmq+2wPU+1mPoDpbX4HofJvz84TMlAMtxA87lzw+Uj85jWIES/PYP5buPYQCgNA9PK764aqWoVoo/UgrHfkxJtxkEriHmhSn67qzCDsg3F/SQCHPxawO50YDLZcpHRZrXIlhDbAcUYpGVYy4eE8YkvyZ1Nfs6QoVxv0yi4Vc/oM6pQ4V3f15E795fS90/my+0s0ygkMPTKNQ7VOht60F9Xb2pc3w/4U3nZGrhgo/BmfDBC84Dkuu6cq/E8xnV5r0SKIZkodtxTnJMPbH8I2XrdRkS+OQ3U/dghpztrH/QJjx90EpPxyRSveQMCNojcEzlE9tu0F+Xn1VmOTFknKjcy6gdcwLv56Xru8LTSVZMOxOpUcZuoWnqfBkygjXE+mf2o16mUOruHSKE5L1PbS2RFHruQ6E/fIiBt5ZSiGuEMPjSexTx+/sSp8JEnOmlNcRcjwtLMP3cjocxWEOM96pKD1fVEHNCSX64oUsNzJkyrSF267LC2YvuQ+7FiYppO88+koxUtMcgXPLrG5J3KQ++V8X4988KK4QfjfjmnMKnyvJhWayZylY0YBFeizI/Ikz2vEQj09pJCUUmzKKKMea3acI73mHUFU5UnwszhUF3V1KX9PHU295PGH20F/WI7Uxh2cOFTpm8VzJXHTcg6xhlJyWtMSN7Ja5PxGKI1XD/U5OOHfhK4AwynDCk9ia7UC/WIod9+cRWEN6prZ9okzQIzJObCySmsvb+i0Id1z2q57slSVyZuillupeR4BF4YYrbxTMMhqsNdsr8gLrBOjC9zKEUXhZGYXDwmRDvaBp8bxX1uTRb6JYVTt1cQ2j0r1OF8benUM+UbjTt1hBhqKc9jU99laIsjwpc3XEnLLyjQOFyXlwF8svTyqEMA8VjVvHLVaX8hlYb9JmUX9jRLMZn1yplNlaKeDia+CBjWYMP4E3fFSk5GBK4Luk/LiopZoOk8eUUSJIGCU7nQR46vMpemLKFHNzrVaZ66lBEoDlFfddT4CxufV3v0KDSwcKE8qk09lZUZd3vzub+UtU3ZPvbQpi1j1TkC9bbetMWKRlxODCHedmxh56xpVQtix+w0bNnfqvMUSlpCbfnVNbiqOu4Rk+twywlzic843Hrye+NswU+xsd1SZ91c44M3bDiZWwOsA2Gz9XekUe1N9gqt7552ZpD9SuXrT1OWZRqxtNNwDOMN52TKMQ9XujOFtDRlzpv6CgYC4dTZ1N/6lM6Vhj9r2k04WEUDT09TOhjC6VpN4fQzJ/DhLHeFjQ10JA+dNYQFlkw2ziEWYVH4aVsXrF0uJSDmAC4IYsvzyusEF+UGqjApnBdUu92DCUBReqSVitFtVL8H6VI2mCkXJgN5hSmpN4tBvrYp5zEcBKAj3DmuMJrFyW4yEH+cpCSpYtXe/MUDiVfi/fMy31EmGKvTcuutKaIjBbCiLT2cJTeofEnw4SBge7UPyOMevuVXv6+Enf4btkoIcQ5krrmTZWKfEz7wlXUxLlYqukwrxUepJdyLJIJn3nx3KfUwOqleolpAicure8/C1/jtCB7EvbTkj2GeeGjU/Rcaiy9dtQmvFpokrOdfJSvzo44gc9l1NnloXrOMqF++jnZy6iX4BPqW32y/f2s3y28EEjBtczU9OgBQfYyDi2jNt5oISTvA6nh3rdwpBB5azKFZfTH/fcR+gXCaFBmTxp7rLcQfrwHDfO1p3HpLYTVN9rRxNQ6tLCglrAyx0Cr4EtwbVJmJ+91ZOieFVOACcSRY/AzipQ8yC0A3+GbQ4p/J/4O16F4m5LK4XhxS4zknWUQ/HPx4iKMSxsVzyaMOzEGqULIcPpeN5zJiw+VDRjPjtyFMlgU+209Qbb/Ti1h8bkGNMlal6L89YTxppfESYrwNxXGZrSmoZY2NMjcXphwrKcUg+/mGihwZd8OvvekRCPTpmQttfungxoX7xOaXMqWAiovWnYInHLwxXQTPWOKFxomJVDDlBRqaHMo3nSqF+eCsBOEl775np7P81PjfxwV+PT3M+Y4evnaT3Jqi+H387kMjrgWUkxy3WdSlBf8JomHeMm0U+DNLd7LeD13p/DmbR+1zFlBob/FCG3NkZI1kK0DExboT5HHe+L+3xSGWNrS2PQ2FJneTBiX2ojGm1+mqb76wiRbHZp/oiHtv/+YwBmQXZhVLE5RiiGXdZgsXCSFFxkz8R4nr1mADOCCXAMbFNsyyHkeZD9HSeETYnFzjZQWYRBcE/GmaLx5ruJcAEGvwoXiFa4v5gloHSrGAevhx4wk8bwS/6WBDnxfFZl14D4coaO1aL6nhhBlqknTLPCafU8I0Rn1aYLteZpR3FQYbG5LRjwtg76cJHCZ517fLKZOZ5YpP+6jZhmr/mtZnLO9NPVvEl53rae/p+2kl1O3CZzW+IUkKIv3gPBC6n56rqCIXrr8rcDphJ63xVbmhuB0CDxzePXq11XH+ErK6DnTvsoQfI64fsm8i162bBP+7uOaYOupiXuD0BRtapS7o2rZunSbLEp1vrhO6Hb2QzKyU1kyVuiF+x1oepNmFjUTIm0vUBSGh+n+J4Vplr9IvwX7cG5uTTpQ/phsjzM7P4cgfzTQ7iNKNmYWsRlV+SccsBw+WP1AimKDi+BeYqD0uYp7Jt4zCfKPVEwLqpWiWin+SCmSZxkpa6BBMIEA/uBiBQG5+JBnNT60V+HipsUYo2z5ykc3dTpU+EBZAKfGdQ/KcEWZBx9jKRo436QssWoBuuUBZa4DzpK/Jk1xPieMcLWmQflhFH4nWujiGal1TFMihe4PKvZKctYJHe74qWn2Bqmmw7QKLKBm1vmSCZ/hJebXXSBzq/CKdy+mnHGSgY55zryfmtwoo6ZXcgXOqsP5IZ5LjaHns33CMykxUuTu9cAWxbOOGrtW0+umxUJTyzxqnbmIOpxYIzTxLZftb97YYloEdC+jje09gdchumWGY7o5TegLH4KHDFYGZnHhozTHWoOWo9+Ypeivpeg3zq/NcH2wufAL9lxUeOhYmPQfOUTgJ5Rg6uk7phw/hYfbqsMG48BDnoWH3TROyYEhSIbcMwcp9jlQipQPYCm6GARTd/gIw/HHSMUPP8O+FH7GLiUJF809apDsJ8xn93QHLvuq4r+qGdgSLis+jG+rs/RAK7MZ3u5WnybqYpZiprIk00DRgTpCRGZLGugJpQn3PhB6ZA6nd68toJ6XFwvtMmZS94cJ9PaVvUIL/0Jqk7dcyisxnTPfo/a2KdTKNFVxz5KCt00TpwuNYU04GQqnJGQ4C13TK3n0Suo64VX7BkkYwqmWGl8/LjT6NIsa+TdT44TpQuv0eXIuo5U5SpAg28Bkan1wrNAmdyk19y6UnU6GFUL2Ms7NEwbeXErdnENkhsEMzOpJY+BDTIN1YGZDIVbCym70KJxsbissdHyBsgQzCi/6NfkfyuJkCP1bAx36l8Lnc3gDLMWrHDunO9/WPYoLD7kP1uHQJCUbimDtBYvRRXFNhVKkRhjpYDODcKwbTIoR1mG0EocPJUOrEjcqFiiFDY06dVFJK9UUSJuSFU6fcxiKEntKWYbG7IRmbzEracV6bG1LkrIJirECT0O09y/CePfrNLa4Fw0oGiCM/X0adTJj+np/pdDjy0XUKj6icllcSjQGZlNn9wShW+pAqabTMXmkwKURQrKnUGvPLKFx6mJ6xb6JGpWlCi8nr8Xws7cyv1T78hxqehzDRMLqylQAL5kxBU5dXnmMr3PGZApJn0hvm0cp8YMkwDbEHyG09kRJcExw65uX8WVjq2LZOsQ5nEb/MpWGYLrJjDnWiyLSmlUOFxw1xUcp1scrPMNz4OnndNjMPvx/KR68fccUtg6Ft6qy9J/6DU//WcihTPHievGYJCRvU1IxdDijMCKMU9L6GehED7ynmeKMqFaKaqX4I6VYMsRIv3UwCIlt0AAMIXH9lWQMIc6FMCmblYBJnZZ/PlRMUJBj3xjo43vKXgh5H24i9yellHTunJqpZOGGCnAj5lxlG661BMPQDPvjwtiU12g475W4Q4XwnyfS8B+r9koG/rKEBvBeyRcLhbbm8VK3s5trqMBF2HrG4/OBfkL/w0Opi2UgdcS0lmmN4abF2QTJYsuwInCeyub564Rmmauow+9pklKozS8BgfNDtPTNp/b2SUKnBM0J0Tezv9DzYBcJw+/uHSywErZJGSfBMRIgAx9i0O3lFOIYLgy+NFmWrsNsoQKvQ4Snvk5T4VAy7HPxkJGcpaTlwY+DcFP8ykYMDyXo18O/KrHo062mqkCoM9cgRwzTPz9UXJCRDw+khR9qYIU8YyHXxMGKrw9eaw/5d1S84VCKxPFGMrc0CDld8TMUwu6nZE+GxcDcNX2j4obGFcE3KChRrpXji/hLixUeyzjN4j6vsg03YCusOnNwmBdQ+AYzlM1wXOcnGmhq4iNChOV5GpXBeyVDhLDUUBr3r8mwGNFCNz51nTaG+n41R5DKvv5xFOYeIHAWnZGONlJeiRl1KJT6+AdQJ8cYobVtKjWxzJW0xkzzjJXUJGUWZjUxQsg/46lp6mx9PVNpdHAmtbJEU2dPhMCHfYdldqMRgY7CKHdbGpXflcK8/YXOnlHSLo6WYroExlAX6yDZ6WR4c0v2Mm4OEWRhyvQyZhi1BHbGuQQDF7JnjvIMLx3+nE3h/twBn2OHQyn7Bb4d+j3njJLi1AL3wSCaowW88wkfcavi+BC/Q677+ihpvWFtMEKYWymp0AeDc5qR7CEGwQdFsL9b5WgmYwgpXgbFWKuk7oRg0wyUn60Ulqo2nvhciYeQub72F78rF/C3RNzQ5e8VWRZHo2NcCnvKs/fBWUpXojPqUGR6c4r6tqfAitEntRsNPjlEmPBwKo17MIPCCkcKXOaZq/r2c3YTBiW/IRX5gvW2+sW1IqO1G/VMGyx0DrxHHU5voJBf44XQh8niuAbTEb59HkNUeYK83uW3BOHtL3bAkXyvMhVAWHKoHPQNHuPjk2/DLC1ogLurwBHXEmCbM1wIvzuDItHu4Z+MEMJYIX4bQjNvhAm8SjnV24DmwUIwPMNgh3KXTZEEI3i4rvykxFg1mPriPcWDh5FzYRacVH56qMoQDJQ6eUQdzKzNSt5iXA/DRtZEJXUAlA4yt3ZR0qOqlaJaKf5IKZJmGqlgmEFwjICwp2CYmKrk4gK+5bj4FqWMzVcMpq1ZyokS3U63uZXij9Dgr2CGbEoKXtuxT2+EsUNhNuL1+CKF/YlNmOLOSXtEmJ5Rm5Z83ZoiAi2EYd52FH6iO0We6iOI8LPCqJert8DF4LvFd6fxZ3sLw52tabSrGQ1xdhC4CFsfvKdnxlCh92dzqJ33fWqXOV3g5Oec67rbkRlC9xNz6E3HBHozI0oSjjGcY6rPpUXUxfyu0Dd7oKQDGB7oJIz2taAhyc3lkA7TeVcIFBDf6ewt9EV7+7neobElvYXxZ3rRUG97CsfnmNU/t5Ol67m5tYRVOXAkXepQMuxD8JBxwKrsjkXfwhE12ZRPvoBcPtZNMMbj0SHj+GEleTeGFzxwnnVKYA1YgP9HK1nsIuDhzxuquHidIneJkWzQFObBKmhOhIHK4xTHXPzEBdO2KHehiV/AF7DvUfhEWSasQ0qCkgwscHoy0hUnlCIeTs7l+8pGzFCOQ5PXeBQvfJIlcDhj7tYSlpxrQJNtdWlaen0h0v6K7JWMczURwjPa0OjMDjTQ9JbQ39xBisFPvzdUmPBRCL13vJOUaGS4Kh8Xt+tq4xJLsDacg+u7RVIagel+aJwkP+/hHqg4QWAk9cgfR119Y4XeF+b+V46piPJoyQ/Bp9+Y6C+70biCjjTjwVBh9NEeci4j2MZR6R0kppKjpZhxuJ/x1ldoWqCBMMlSh+Ydb0ix5Y8Jc926DsFrOAw7lew7fF6umGFtOWkd72cwFrzHnFglA97A9EJxbLuU83BO75+CLDcq5fhs8my8Fqs438PPlXi4hyllnMU/D0qR2tcg/I7pig8m5aetysnVMFEQ/HWnwkcLbVCcou0Ka54DJsmyQ3GgEW5YhuTtCh87dMP8nbuhmKDBvq/R+M+UfZ9Aaa5WLYvvv/coLT72KM3B08FEW2rRdPdjNM37pPC+pwFNsL4ARWFleYWiS5pJctEJp7sIIxxa6nnOv/sJw0t7S0W+YL0tFmqng/2pf+5ggUsjcBLYgb6ulfR3v0N903qBd4Wu5gEy8wmmE+JdXN7d/bB8oDAosTmNy2lD44o7CcNgBfiQTjAEnyOu3/c2pGm+pwSOSOPgmODW94LDNXVjK1eXrXdfNMii1N5jCs8wdqAvLt5VUmyqDPb9ihc4durQztghD/d6WPzNihf/t0Nu39uVq3hoC5Ya6PZexRQOY7AC1+2rlPGGWLVSVCvF/1WKRUb6YZJBODLeQL/Ah0ifptyFQphmwcFJUh6gMUe3YdyZqQQw5zXPx00sUnhJ/CCPV7sUGxpwHmMenz1gzt/TvZKCB8piTFNdGFb2fa0sgmlbBl9lmUNZmKJ+BwenMvM8NWh68iM0y/uoMMX2JEWXtaSo822E1Q870uJf36ZRWW8KM+4NoV6mrlKRj4l4GE1Dz4dLuiBmsLkNjfG2pMjsNyoZ62tJQ1Jb08jcUIHXJDgt4fibk4Seid0o6vYQGlfaSZjzTQdaTW/RzAttBT7bOcFURw7pMNGJNSUMf3mmQeCYSkmJbFFWsjKgX3d9ofBeBoc05v+uHDis085gXTfuTx4yPPsVCyYCLsgsaZFimgfZQS7+2UrWWh0yrpqVlBkG+jf8wqIPle/gN378PhinlC2FUmTPN9I9vMiY4InexZfkva9cwgXvYtw5ukKJj0IjpqMRkxUnZigJE/D/aIXX1JPRoOAG2gHgxywlzqx8gTExBh2Qe13hvZKF6IDkK4oTCrLhkM5KGEkMH6jaB1gExeQzDhvzlVmuWjTe8hwt+6G5sO7uG5I4ZbSnpRB1tTdF/3uwlGhkekCgfVK70hBrW4Er6kQk16Mob20hOq02RabUoyn+12iopaXQ396VeqRwjqkeAicd4xxTwXwRXHx+/f2WtOzH5kKE/UWa7Xu8spDbIoz3fFiHE4cwOyCgXXgwgoV0VqapIph+VLg/eC+DVymZLSm6DsEHfJgTcCrNB9UqM6wQ3Pfu6UoSHm57hMqGcUMx4vCAH16m3IZC/AAL4puk3IdcLYMg5w+U44uhFKvGG+l7mBDmawg4bywUZLmSCSckF5r1y3alPB5P+DzdSWXijbAGww1U+r7igtIkQSk42xpjgxA9cHTOfqW4MD36+H7Vsjiveh6C9dh/XFmCDuD6IbtcCk+3JEo8SdkGxdoPJVnnVhaClWefoPllDYXIFEztPmlUmW5plOsNGn+sU2XNTi7ROIGr8nkbC9Ew59HJf6I5EBrzIdrEmfA58flE/+sCp2eS8o4VKQknYpjihGPBdEILLzWm8YkNaO6xhsLK80/JQd+VDoVPbW2HZdieqPD98LR9n0dhJeHgmODWNy9b88bWmniFVyl5USo45bTinpPjq4YLtg6sDEVTFDvkkdAP1iKiAvytHA7ozZ0KWwcPXIQ7y5UT+P85yP7mBGXLhGqlqFaKP1KKnNlGMnc1CPdw0aMjYWbGKawYP2DMyZismKAIHGdhH6R4+0Bx+mOcGqzE4zN2KI1ro5LFi1nooBvlCh8cOv4NFOSOEoMO2o/hJfOqUszDC4aFRLfCOTIcXlzbrsSh83hzKFhvZIWT901q0IGHfxH20d9oRm4DmnXiNYHTOc75ur1U9WWG2lrTCFsLOH/PC9PNNSW+cZ1V2ejUeMeZtj/TBEsDYYyruaQ2ZieWmXGhvRSZ54RjzBRXA3xvXYp5+IQwy/8nSQfAFZQYPtsZi3bzIR2Gz2VwGH4wiQjHVHIIXc5PCocY7IYynbqlnMEwy0vXwb0MXofgaadts2KaXzFkjFCOQA72nnjfQCV9gsotbbLy7UIdMo5PUHLxmfsYNmzvKEUfsqM5w0jfQrCMHxe7C8X4MlKx4KJncKGbi5UHmI0cxoV8/ZX0btBWXCjBqHjZgqCRvjWKEz5FSRE+c1ThvZIU3HDmSeVzjKXHb0AB8hQucxQPofCiF2OCZTgKJeATT0wAT1UMOnifT+GT7htx/dkOZSWuuf/hU7TkzNPCJNezUgyeyzwzq8s70MzTLWim76/CEr+ertqQrKyL1QIqXC9jTtrjwtSjjWnVg7doXXl7gROXvud/vrLAPCcdW/dpDZrt5MPVBlqdrxYiuDfBez7ZuVX3cDgHyoGfqW4lGfe0Hj5CbKZSck33kHLPKry59b/3MngdwrNeYafSNl6tA2PrqXJxGZXAGFgKyO32UoWdSsfgKsfyPn7PwsN9Lkwpm80bYpOM9O93DcIlXNQOIf8YrtyDBn2GL7T0UzJhmnKH4OZCFVcf1UZfuGLFl2TDUhRsU8wYcjLR6XlZStFxPXp47BMl0ake9blfFQ40tUI57G7lNBTpHDrh+gUlAWYz45BaC2Y3hLcKnWovV7ac0antvlu1hD33/yaCm2RrIEzLfZGmeOrTdMufhWWwPgcgpNxPlWPfG6SaDlc7mpH6JyHaX08+N9nZUJh9+Gnafbcu7f/9UWEOW5rSqjZswv3shoWz5CputC8WQrx6QbmIaeaXp9B+p5Lo0/v+5JbC0VIpsAYFZQpbh/9ctj5eoItSwSknB9s6IquGCw8suLsvvrOLUjgUMoE1sPRXTo9Vh/JGhOLqgQcfynBvgOKZUq0U1UrxR0rhnGgke0uDEHgLSgBBB7op/u4wVb3w5e8qHw0z0O9wPkvZuQQe+B8ZPLWZpmRAIdJWoCE7lFKYfzcEWZKnnECHuSFsD6/Fg+IzesOugMLHB7KhCElmpRDjMZ+FDNZiz+BhBj4Gl7RiuB4J5+Ba7lRW4/eYC5oOgZnt4WJ3NWlB8eMCF4Pn2t8L8DrDhVy5bmew3tZ2XCOuREssLctWNv6/WrT87P9I4nNmVkZNSW28yK/EfKFpCYNrK5ug/AcKcT8lSjKuxcHOaZmKF9/DZTSCCdNZ8HzfHGDLcEwlh9BxcAzD8RA8ZJwsUngv43yaLkox2esMkkI7MENJw/AegOCPDFXYTzwHpTg8QPH1Bt11iGE+h7yz3kY/t1Q87GgmjDSS91WDkNQIf2iKC7eqoD2edrC/nZKLixTAMgQdTd4rMcOSPIhVTGjUPWi+d5PCeyVf40lM3aPwXokfY+zBGMUKwTvgjHpZ0CATnZB/VB0xhhdpzn9W5WRZ8J5P70KQeUruj3gNHbjdoaxDJ29CJ24OKKu9FQGvCcpyCHMxXucyz4ztVzh6BVqNj2GHj4uwcc0trqbDzPfoIad5ScqcOIPkut6coWx06fkXVijmIJSs8Cae0iyFC/SxX8EWkvkFs4nt2w2VxfyOHNcVSg5eYjjI1o73x8coHCDjw4Nl3avw5hbvZdjXK7wwFTdVZxiMCU9/+Sq1Dsw5PLgBthrtlRTgexMK0VpxNoPcX8drrymWkbxLOtBI1joGwdwADX0WAnlZsbyChuGNgRZKIhQjHZpl769c46ECpujqZuXkGgNdRqddcygObjicnLytCgd5cPRP0haFFYWPDqQnK5yu8dgRCP6UcuefBjk2f/W2cvRjg6RbOvu7Eg+Fi0MnmQMKd76cvs5UNuKza2I1aIXZDSGthTIF637zyiJX9bWcUXj1MOemvm7/Tdl4RBfQtjqVNXGa+JzD6JmD6fq9B91Kcq4GM5/8VTmOIYlrvf50X2Hn+cdvqo7x8aktJxQ3GIKfmaBb3cHwOe4v7rfsTQqHMDg2/ceytVkXpW7uUtImqUOZ2l9h63AASpDxhhIHeZpgADJfVvzP4Tsg99S6inNgtVJUK8UfKUWi0UjmPxkEW02MYY9ijPofxQ5FOYIPpbyquHBxe3cDfTNROYbx6yc0IIP9CXAPAjDPxN8OKmzaSnBjmfMU5xyDxAhaligumDkncK1V/Bhe+Hg8H2oWruPGL6FTDymc4i8B5vXoJYWnbuyoZh9SrBC8Febc5Fa4xioXTUvwK1xcdy8Ez7klGT+uz7W/uYArsx1/D9Z238hDA+CaW1xiiaeODOf5SsB1OT6V4VzXNijLoULFlKXtOnFF4ff/gCGv9JTyOYY7zg/BB32ZPLc+GHxIh+FzGdwvHGAroL84hC5vkRJYXbGXYVKSMWTfPlC1l8HrEDztPDNWYR+Ch4zY15QCyNMGubqfVFyP4fq18D1/VszQB8MB/LPPYBBcwA+8FTihKM7auPArSiYsRaC7LnYw5nfxNEMpcqYovFdyB4pRvELhqgAuKAIvnjDpk6Gp+GmdpvhnqrLwghfDjlMpOinfo5R/r4qRmaYEc3DlHVd4L4AFY7UqJfn4Dh6XHUosnur4BA0YZthJ3QOhFT1QHJ9qMXgu88ychoJsxTUtEN6RewrvwSTjmlxNh+ECKgdwXZtNCaBdR6A4XIeL4eRjvM+TWaywQnwNJbClKuXXVRmOZyuHcQ3JdrxY4YM6vIfEEdcMB9lyTGXyWIUDZFJmw29ZprBC/AgrzauUDC9M2QdX7WX4eqgPwdaBscMa+B6BQhiUNOABeytgI2HYiH824RcmAaSAXRWk8wf/ip8vKq4W0Go4m9+OUf73sjgrBi+LB03ZA3QeL4vnTFL2QGsTBuEa4UruRJ3G5sxX+MDrCQj3cqGSuBudiifr3jWFS00UoKOv31HY+Sw8htc9SjqeuACcwfidSgI+nw/znpSk8Ipqam7VMnoWZg6fQMlO/KCUfK1VffMu4yHxK4kYOpKcqgxMIYYsLp4SLI2Qnoj7gNkPFmbzZWhI3I37SulH+E4o3p2rynUMf3xY+5Ns5QKUlA9zB4/x8Vb2IShBMASfI645wDZrguKOUqc+uPXN1sE7qWrZmq03L0oFp5w8w0hrXeUSiIWAXHdWwPI+CLZVsLVaKaqV4g+VYjX+WY5fmLVgCTBX4ANu+BiOZxSerp7qCTMYqtx57//ulVxbqsuuDB9c5aCdpP7KIeALM8gBFKZwOH7ic+YpSiqbzq0wr6eUe2fR0TD/54sVTuf4HQQZTPfn9ugwkRqrpMCMWtfAzFbghNA4GXmGV3EBE4aAy6QUfQVFgdK4spRAgUEK9h763CBFXBnOE8qBQpyznPFD0T1QOPM6JRVjvGk9vitO4eRxXrw/JUVhZ5LTEl76SHGinb+fqbrHdLTRugR9NV3JjsQDMriqj/hcBofhBx1HPnbBIXSsPMwVDDn30e9lE5Xc4bp07e6h8DoETzt9zynsQ/CQwcrAzAPrwMoK1rFSLMI/7+EXhl9kxYitQJQD44/paSX57/giOC1f9VP8vXRx5GKkwgdVT07WcU1YiScbmps5UMnqiuvBLykNVZKgxakDq8bDNFiVUmi/bb1yBQIsP40n2K3YD8LxRCdzZmDm/lVNxpaWoHABm/3oEOuHyuEtKrCkvQpvzgXy1flkeFbwPWYhv5RXwT5AAYQVa1bS89RXSNyrpOCah+E8W+cq+3H/OVBk7wHlcLq2K9jGny+gz2Cl8q3KAyjENbTBvlEpwTV9uG/HRIX70IIH71h3hQ/q8LkMjrRn3KO0bzlaimGnMnWw7nQyrBDZUISzYUr227oOYW6osFPJPgRbB4YVgg3B+xUsYaWYj38W4BdmAlgMVlcQAxw11NkUh/NlaNobuiTOXEYjzWjw9fEKL4tfiIQg+iquIXhCQWKIktLGQL+9hb81VzKgJHzANSNScWPoCeBG72UqOXiKsjAE3DtZAYaSfJhi634ll520GHTuJiUJliYfJIQrnISFTbN7j+LF+z1pVSWcv/1ZnUezWeEnm3dmv/xH1Xu4kgGXuchMVDL4u2YZZJGIyYPyJcOZDi4meWEtcqFEKbsU/szvsHgPTimHYSl4Ye9OpnIE1s2PawRTAWSgP/igr6W5crcjvqclHpIQJX8o2jOsqo8/DleH8vIohWeH56AUwWVrRyuDLEpZ6yo8w9hjqBou2BCwMrDFYNZUK0W1UvyhUmxetoxW4D/MQsA+xvYKZA0DpA6oYISRMicZyVpB2Swjlc410sYIpXSRUWL8OPhTWGGkfLzmnqqkjDPSrpFG8kcqtvdx/ZlGip+jJOKzjk1GOmmpwGqktDh8body0m+U3NFlGcqJHCMdw8+kzYp9o5FiFxqpcL3iWGekZFxv52qFSzbHxRpp1SplBdpXeBjXKq7i6BEjrcPnuBA8s2+vfi5hh8JFc7y4ZmClkoDvs+K7nLsUznVekov7dijcXm63fV8Fe/W+OOEYwzmmkpeB+YoN/eBEX2VPVGLQXxb0W1qU4pyj/VrZx+gzDqFbH6FwPARvf/PGljABfYxrOAYpLM8EsKUC9iF4yFhZwV7ow/8HkQijoXc9bVkAAAAASUVORK5CYII=>

[image4]: <data:image/png;base64,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>

[image5]: <data:image/png;base64,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>

[image6]: <data:image/png;base64,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>

[image7]: <data:image/png;base64,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>

[image8]: <data:image/png;base64,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>