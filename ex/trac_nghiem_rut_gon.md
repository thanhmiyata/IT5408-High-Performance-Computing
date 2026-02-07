# Đáp án Đề thi IT5408 - High Performance Computing

## I. Multiple Choice Questions

**Q1. Which of the following is not type of Flynn's Classification of computer?**
- A. SISD
- B. SIDD
- C. MISD
- D. MIMD

**Answer: B. SIDD**

**Q2. In shared Memory**
- A. Changes in a memory location effected by one processor do not affect all other processors.
- B. Changes in a memory location effected by one processor are visible to all other processors.
- C. Changes in a memory location effected by one processor are randomly visible to all other processors.
- D. None of these

**Answer: B. Changes in a memory location effected by one processor are visible to all other processors.**

**Q3. In distributed Memory**
- A. A computer architecture where all processors have direct access to common physical memory
- B. It refers to network-based memory access for physical memory that is not common
- C. Parallel tasks typically need to exchange data. There are several ways this can be accomplished, such as through, a shared memory bus or over a network, however the actual event of data exchange is commonly referred to as communications regardless of the method employ
- D. None of these

**Answer: C. Parallel tasks typically need to exchange data...**

**Q4. Synchronous communication operations referred to?**
- A. Involves only those tasks executing a communication operation
- B. It exists between program statements when the order of statement execution affects the results of the program.
- C. It refers to the practice of distributing work among tasks so that all tasks are kept busy all of the time. It can be considered as minimization of task idle time.
- D. None of these

**Answer: A. Involves only those tasks executing a communication operation**

**Q5. Granularity is?**
- A. In parallel computing, it is a qualitative measure of the ratio of computation to communication
- B. Here relatively small amounts of computational work are done between communication events
- C. Relatively large amounts of computational work are done between communication / synchronization events
- D. None of these

**Answer: A. In parallel computing, it is a qualitative measure of the ratio of computation to communication**

**Q6. Parallel Overhead is?**
- A. Observed speedup of a code which has been parallelized, defined as: wall-clock time of serial execution and wall-clock time of parallel execution
- B. The amount of time required to coordinate parallel tasks. It includes factors such as: Task start-up time, Synchronizations, Data communications.
- C. Refers to the hardware that comprises a given parallel system - having many processors
- D. None of these

**Answer: B. The amount of time required to coordinate parallel tasks...**

**Q7. Scalability refers to a parallel system's (hardware and/or software) ability**
- A. To demonstrate a proportionate increase in parallel speedup with the removal of some processors
- B. To demonstrate a proportionate increase in parallel speedup with the addition of more processors
- C. To demonstrate a proportionate decrease in parallel speedup with the addition of more processors
- D. None of these

**Answer: B. To demonstrate a proportionate increase in parallel speedup with the addition of more processors**

**Q8. Non-Uniform Memory Access (NUMA) is?**
- A. Here all processors have equal access and access times to memory
- B. Here if one processor updates a location in shared memory, all the other processors know about the update.
- C. Here one SMP can directly access memory of another SMP and not all processors have equal access time to all memories
- D. None of these

**Answer: C. Here one SMP can directly access memory of another SMP and not all processors have equal access time to all memories**

**Q9. Task dependency graph is?**
- A. directed
- B. undirected
- C. acyclic directed
- D. acyclic undirected

**Answer: C. acyclic directed**

**Q10. In the threads model of parallel programming**
- A. A single process can have multiple, concurrent execution paths.
- B. A single process can have single, concurrent execution paths.
- C. A multiple process can all tasks are kept busy all of the time. It can be considered as minimization of task idle time.
- D. None of these

**Answer: A. A single process can have multiple, concurrent execution paths.**

**Q11. Data dependence is?**
- A. Involves only those tasks executing a communication operation
- B. It exists between program statements when the order of statement execution affects the results of the program.
- C. It refers to the practice of distributing work among tasks so that all tasks are kept busy all of the time. It can be considered as minimization of task idle time.
- D. None of these

**Answer: B. It exists between program statements when the order of statement execution affects the results of the program.**

**Q12. Which of following is an example of domain decomposition?**
- A. Matrix multiplication
- B. Merge sort
- C. Quick sort
- D. 15 puzzle

**Answer: A. Matrix multiplication**

**Q13. In designing a parallel program, one has to break the problem into discreet chunks of work that can be distributed to multiple tasks. This is known as**
- A. Decomposition
- B. Splitting
- C. Compounding
- D. Computing

**Answer: A. Decomposition**

**Q14. Functional Decomposition is?**
- A. Partitioning in that the data associated with a problem is decomposed. Each parallel task then works on a portion of the data.
- B. Partitioning in that the focus is on the computation that is to be performed rather than on the data manipulated by the computation. The problem is decomposed according to the work that must be done. Each task then performs a portion of the overall work.
- C. It is the time it takes to send a minimal (0 byte) message from point A to point B.
- D. None of these

**Answer: B. Partitioning in that the focus is on the computation...**

**Q15. Which value has the speedup of a parallel program that achieves an efficiency of 75% on 32 processors?**
- A. 18
- B. 24
- C. 16

**Answer: B. 24**

**Q16. The type of parallelism supported by GPU CUDA is best described as**
- A. MISD (Multiple instruction single data)
- B. MIMT (Multiple instruction multiple thread)
- C. SISD (Single instruction single data)
- D. MIMD (Multiple instruction multiple data)

**Answer: D. MIMD**

**Q17. Within a parallel region, declared variables are by default**
- A. Private
- B. Local
- C. Public
- D. Shared

**Answer: D. Shared**

**Q18. Collective communication**
- A. It involves data sharing between more than two tasks, which are often specified as being members in a common group, or collective.
- B. It involves two tasks with one task acting as the sender/producer of data, and the other acting as the receiver/consumer.
- C. It allows tasks to transfer data independently from one another.
- D. None of these

**Answer: A. It involves data sharing between more than two tasks...**

**Q19. ________ is an object that holds information about the received message, including, for example, it's actually count**
- A. buff
- B. count
- C. tag
- D. status

**Answer: D. status**

**Q20. In a CUDA program, calling a kernel is typically referred to as?**
- A. Kernel invocation
- B. Thread invocation
- C. Block invocation
- D. Grid invocation

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
- A. 1000
- B. 100
- C. 10
- D. 10

**Answer: A. 1000**

**Q22. If both export OMP_NUM_THREADS and omp_set_num_threads both are used in a program than the threads will be set to launch the parallel section in accordance with?**
- A. OMP_NUM_THREADS
- B. Whichever sets minimum number of threads
- C. Whichever sets maximum number of threads
- D. Omp_set_num_threads

**Answer: D. Omp_set_num_threads**

**Q23. Which of the following correctly describes a GPU kernel?**
- A. A kernel is part of the GPU's internal micro-operating system, allowing it to act as in independent host
- B. A kernel may contain a mix of host and GPU code
- C. All thread blocks involved in the same computation use the same kernel
- D. All of the options are correct

**Answer: C. All thread blocks involved in the same computation use the same kernel**

---

## II. Fill code

**Exercise 1:**
- 
**Answer:** `*d_res = (f_device(x0 + h) - f_device(x0 - h)) / (2.0f * h);`

**Exercise 2:**
- 
**Answer:** `for (i = index; i < N; i += blockDim.x * gridDim.x)`
