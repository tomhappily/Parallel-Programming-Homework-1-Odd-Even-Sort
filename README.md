# 2024 Parallel Programming HW1 Report

## Implementation

### Handling an Arbitrary Number of Input Items and Processes

1. **Data Partitioning**:

   * Calculating each process’s portion of the data by dividing the total number of items (`global_count`) by the number of processes (`global_size`).

   * Using `global_base_size` as the base size per process and `global_remainder` as the remainder, we assign one extra element to the first `global_remainder` processes. This ensures a nearly even distribution, even if `global_count` isn’t perfectly divisible by `global_size`.

   * Each process then calculates its offset position (`local_offset`) to know where to read and write its portion from/to the file.

     ```C++
         global_base_size = (global_size > global_count) ? 1 : global_count / global_size;
         global_remainder = (global_size > global_count) ? 0: global_count % global_size;
         int local_size = global_base_size + (rank < global_remainder ? 1 : 0);
         int local_offset = (rank >= global_count) ? 0 : rank * global_base_size + min(rank, 			 						global_remainder);
     ```

2. **Process Grouping**:

   * Using `MPI_Comm_split`, the program creates a sub-communicator (`g_active_comm`) with only the required processes for handling `g_count`. This avoids having idle processes if `g_size` exceeds `g_count`.

     ```C++
     int color = (rank < g_count) ? 0 : MPI_UNDEFINED;
     MPI_Comm_split(MPI_COMM_WORLD, color, rank, &g_active_comm);
     ```

     

### Sorting Methodology

1. **Initial Local Sort**:

   - Each process reads its assigned data from the input file and performs a local sort on its data using `std::sort`. This initial sort ensures that each process starts with an internally ordered subset of the data.

     ```C++
     if(local_size != 1)
     	sort(local_data.begin(), local_data.end());
     ```

2. **Odd-Even Sort Phases**(the sorting algorithm proceeds in rounds of *odd* and *even* phases):

   * Valid processes within the g_active_comm execute `odd_even_sort_process()`.

     ```C++
     void odd_even_sort_process(std::vector<float>& local_data, int rank, int local_size);
     ```

   * During an odd or even phase, processes with the corresponding rank compare and exchange data with their neighboring process, as long as the neighboring process exists.
   * The signal is set to true if the process needs to swap data with its neighbor, and this signal is sent regardless of its current value.

3. **Merging Neighboring Data**:

   * When two neighboring processes detect that an exchange is needed, they use the `mergeSortedArrays` function. This function merges the data from both processes into a buffer, ensuring that the merged data remains sorted.

     ```C++
     void mergeSortedArrays(std::vector<float>& arr1, std::vector<float>& arr2, std::vector<float>& buffer);
     ```

   * After merging, each process retains its half of the merged data to maintain balanced distribution across processes.

4. **Checking for Completion**:

   *  `MPI_Reduce()` reduces the `signal` variable from all processes in `g_active_comm` to a single value stored in `buff` on the root process (rank 0). A *true* `buff` means that additional sorting rounds are necessary.

   * `MPI_Bcast` broadcasts the value of `buff` (the result of the reduction) to all processes in `g_active_comm`.

     ```C++
     MPI_Reduce(&signal, &buff, 1, MPI_C_BOOL, MPI_LOR, 0, g_active_comm);
     MPI_Bcast(&buff, 1, MPI_C_BOOL, 0, g_active_comm);
     ```



## Experiment & Analysis

### Methodology

1. **Performance Metrics**：

   * **Process**:  The number of processes used in the experiment.

   * **Total Time**: The entire execution time from reading the data to writing it back after sorting.

   * **I/O Time**: The time spent reading the input data from the file and writing the sorted data back to the file. This is measured using timestamps before and after the I/O operations.

   * **Communication Time**: The time taken for processes to exchange data during the sorting. This includes the time for sending and receiving messages as well as any synchronization.

   * **CPU Time**: This metric represents the time spent in actual computations and sorting processes, excluding I/O and communication times. 

     

### Plots: Speedup Factor & Profile

1. **Experimental Method**:

   * **Test Case Description**: The size of the dataset chosen for this experiment is **536,870,888** elements. Due to the moderate execution time of this dataset, it allows for more accurate measurements. 
   * **Parallel Configurations:**: Since the sort function is used to sort the data within each process before the odd-even sort, the chosen numbers of processes (4, 12, 24, and 48) provide a better demonstration of the impact of the number of processes on speed. (Using a single process would require the entire dataset to be sorted with the odd-even sort, resulting in a **significantly longer** execution time.)

2. **Analysis of Results**:

   ![time](https://raw.githubusercontent.com/tomhappily/Parallel-Programming-Homework-1-Odd-Even-Sort/main/img/Figure_2.png)

   ​												 	 **Figure 1**

   

   ![time](https://raw.githubusercontent.com/tomhappily/Parallel-Programming-Homework-1-Odd-Even-Sort/main/img/Figure_1.png)

   ​													  **Figure 2**

   

   * **Figure 1**: It shows the performance of different numbers of processes in terms of total execution time and the speedup ratio with 4 processes. The results indicate that as the number of processes increases, the execution time does decrease; however, the reduction is not significant and does not meet the expected results.

   * **Figure 2**: It shows that as the number of processes increases, the speedup of CPU time has nearly reached the ideal expectations. However, the I/O time remains almost unchanged, while the communication time has increased.

     

### Discussion

1. **Bottlenecks**:
   * **I/O performance**: The primary bottleneck in the system appears to be I/O performance, primarily due to the fixed nature of disk read/write times. As the dataset size increases, the time taken for file operations does not decrease significantly, leading to longer overall execution times. 
   * **communication time**: Communication time becomes a concern as the number of processes increases. The overhead associated with synchronizing and exchanging data between processes can further contribute to the overall latency, particularly in scenarios where many processes attempt to communicate simultaneously. 

2. **Improvements**:

   *  **For I/O**: To enhance I/O performance, consider using faster storage solutions such as SSDs, or implementing more efficient I/O strategies, such as asynchronous I/O or parallel file writing techniques that can allow multiple processes to write concurrently.

   - **For Network**: Optimize the communication patterns to minimize the amount of data transferred between processes. Techniques like message aggregation and reducing the frequency of communication can help alleviate the impact of network congestion.

   
## Experiences / Conclusion

### Conclusion of the Assignment

 This assignment provided valuable insights into the implementation of parallel sorting algorithms using MPI, specifically the odd-even sorting method. By analyzing the performance metrics, we were able to identify bottlenecks in I/O and communication times, which informed our understanding of the complexities involved in parallel processing.

### What I Learned

 Through this assignment, I gained practical experience in parallel programming and the challenges of scaling applications. I learned how to effectively measure various performance metrics, including total execution time, I/O time, communication time, and CPU time. Additionally, I developed a deeper understanding of how different factors, such as the number of processes and data size, impact overall performance.

### Difficulties Encountered 

1. **Debugging Parallel Code:** Debugging MPI programs can be more complex than serial code, as it requires synchronization between processes and careful management of data communication.

2. **Performance Measurement:** Accurately measuring and interpreting different performance metrics, especially in the context of concurrent execution, posed challenges. Ensuring that all timing measurements were taken correctly and accounted for was critical but sometimes tricky.










