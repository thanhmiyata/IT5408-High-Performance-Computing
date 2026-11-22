# Đáp án Đề thi IT5408 - High Performance Computing

## I. Multiple Choice Questions

**Q1. Which of the following is not type of Flynn's Classification of computer?**
**Answer: B. SIDD**

**Q2. In shared Memory**
**Answer: B. Changes in a memory location effected by one processor are visible to all other processors.**

**Q3. In distributed Memory**
**Answer: C. Parallel tasks typically need to exchange data...**

**Q4. Synchronous communication operations referred to?**
**Answer: A. Involves only those tasks executing a communication operation**

**Q5. Granularity is?**
**Answer: A. In parallel computing, it is a qualitative measure of the ratio of computation to communication**

**Q6. Parallel Overhead is?**
**Answer: B. The amount of time required to coordinate parallel tasks...**

**Q7. Scalability refers to a parallel system's (hardware and/or software) ability**
**Answer: B. To demonstrate a proportionate increase in parallel speedup with the addition of more processors**

**Q8. Non-Uniform Memory Access (NUMA) is?**
**Answer: C. Here one SMP can directly access memory of another SMP and not all processors have equal access time to all memories**

**Q9. Task dependency graph is?**
**Answer: C. acyclic directed**

**Q10. In the threads model of parallel programming**
**Answer: A. A single process can have multiple, concurrent execution paths.**

**Q11. Data dependence is?**
**Answer: B. It exists between program statements when the order of statement execution affects the results of the program.**

**Q12. Which of following is an example of domain decomposition?**
**Answer: A. Matrix multiplication**

**Q13. In designing a parallel program, one has to break the problem into discreet chunks of work that can be distributed to multiple tasks. This is known as**
**Answer: A. Decomposition**

**Q14. Functional Decomposition is?**
**Answer: B. Partitioning in that the focus is on the computation...**

**Q15. Which value has the speedup of a parallel program that achieves an efficiency of 75% on 32 processors?**
**Answer: B. 24**

**Q16. The type of parallelism supported by GPU CUDA is best described as**
**Answer: D. MIMD**

**Q17. Within a parallel region, declared variables are by default**
**Answer: D. Shared**

**Q18. Collective communication**
**Answer: A. It involves data sharing between more than two tasks...**

**Q19. ________ is an object that holds information about the received message, including, for example, it's actually count**
**Answer: D. status**

**Q20. In a CUDA program, calling a kernel is typically referred to as?**
**Answer: A. Kernel invocation**

**Q21. Consider the following piece of code:**
```c
int main(void){
    omp_set_num_threads (10);
    #pragma omp parallel {
        for (i=0; i<100; i++)
            printf("Parallel Processing\n");
    }
    return 0;
}
```
**How many times "Parallel Processing" will get printed?**
**Answer: A. 1000**

**Q22. If both export OMP_NUM_THREADS and omp_set_num_threads both are used in a program than the threads will be set to launch the parallel section in accordance with?**
**Answer: D. Omp_set_num_threads**

**Q23. Which of the following correctly describes a GPU kernel?**
**Answer: C. All thread blocks involved in the same computation use the same kernel**

---

## II. Fill code

**Exercise 1:**
**Answer:** `*d_res = (f_device(x0 + h) - f_device(x0 - h)) / (2.0f * h);`

**Exercise 2:**
**Answer:** `for (i = index; i < N; i += blockDim.x * gridDim.x)`
