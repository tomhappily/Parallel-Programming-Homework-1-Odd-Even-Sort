#include "mpi.h"
#include <vector>
#include <string>
#include <algorithm>

// #include <time.h>
// #include <cstdio>
// #include <numeric>

#define Single_Data 0
#define Bool 1
#define Multi_Data 2

int g_count, g_remainder, g_base_size;
int g_size;
char* g_in_file, *g_out_file;
MPI_Comm g_active_comm;

//total time
// struct timespec start, end, temp;
// double total_time_used;
// //io time
// struct timespec i_start, i_end, i_temp, o_start, o_end, o_temp;
// double io_time_used;
// //comm time
// struct timespec c_start0, c_end0, c_temp0, c_start1, c_end1, c_temp1;
// struct timespec c_start2, c_end2, c_temp2, c_start3, c_end3, c_temp3;
// struct timespec c_start4, c_end4, c_temp4, c_start5, c_end5, c_temp5;
// double c_time_used;
// std::vector<double> comm_time;


void odd_even_sort_process(std::vector<float>& local_data, int rank, int local_size);
void mergeSortedArrays(std::vector<float>& arr1, std::vector<float>& arr2, std::vector<float>& buffer);

// void temptimeCalculate(const struct timespec & start, const struct timespec & end, struct timespec & temp);
// void timeCalculate();

int main(int argc, char **argv)
{   

    // clock_gettime(CLOCK_MONOTONIC, &start);

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);

    g_count = std::stoi(argv[1]);
    g_in_file = argv[2];
    g_out_file = argv[3];

    //避免多余process
    int color = (rank < g_count) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &g_active_comm);

    if (g_active_comm != MPI_COMM_NULL) 
    {   
        if(g_size > g_count)
            g_size = g_count;

        g_base_size =  g_count / g_size;
        g_remainder =  g_count % g_size;
        int local_size = g_base_size + (rank < g_remainder ? 1 : 0);
        int local_offset = rank * g_base_size + std::min(rank, g_remainder);
        std::vector<float> local_data(local_size);

        // clock_gettime(CLOCK_MONOTONIC, &i_start);
        MPI_File input_file;
        MPI_File_open(g_active_comm, g_in_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file);
        MPI_File_read_at(input_file, local_offset * sizeof(float), local_data.data(), local_size, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_close(&input_file);
        // clock_gettime(CLOCK_MONOTONIC, &i_end);

        //先对每个process的data进行排序
        if(local_size != 1)
            sort(local_data.begin(), local_data.end());
        
        odd_even_sort_process(local_data, rank, local_size);

        // clock_gettime(CLOCK_MONOTONIC, &o_start);
        MPI_File output_file;
        MPI_File_open(g_active_comm, g_out_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
        MPI_File_write_at(output_file, local_offset * sizeof(float), local_data.data(), local_size, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_close(&output_file);
        // clock_gettime(CLOCK_MONOTONIC, &o_end);

        // clock_gettime(CLOCK_MONOTONIC, &end);
        // timeCalculate();

        MPI_Comm_free(&g_active_comm);
    }

    // std::string filename = "./time/" + std::to_string(g_size) +"process/rank" +
    //                      std::to_string(rank) + ".json";
    // FILE *fp = fopen(filename.c_str(), "w");
    // fprintf(fp, "{\n");
    // fprintf(fp, "  \"total_time\": %f,\n", total_time_used);
    // fprintf(fp, "  \"io_time\": %f,\n", io_time_used);
    // fprintf(fp, "  \"comm_time\": %f\n", c_time_used);
    // fprintf(fp, "}\n");
    // fclose(fp);

    MPI_Finalize();

    return 0;
}

void odd_even_sort_process(std::vector<float>& local_data, int rank, int local_size) 
{
    float rightmost_neigh;
    float leftmost_local, rightmost_local;
    //rank 0中的signal控制sort是否结束
    bool signal = true;
    bool buff = true;
    //count为0时为Odd Sort, 为1时为Even Sort
    bool count = 0;
    int first = 0;

    int next_size = rank == (g_remainder - 1)  ? local_size -1 : local_size ;
    std::vector<float> buffer(next_size + local_size);
    std::vector<float> recfrnext(next_size);

    while(buff || first == 1)
    {   
        first++;
        signal = false;
        buff = false;
        leftmost_local = local_data[0];
        rightmost_local = local_data[local_size-1];
        if( ((rank+count)%2 == 0) && (rank + 1) != g_size)
        {   
            // clock_gettime(CLOCK_MONOTONIC, &c_start0);
            MPI_Recv(&rightmost_neigh, 1, MPI_FLOAT, rank+1, Single_Data,
                     g_active_comm, MPI_STATUS_IGNORE); 
            if(rightmost_local > rightmost_neigh)
                signal = true;
            MPI_Send(&signal, 1, MPI_C_BOOL, rank+1, Bool, g_active_comm);
            // clock_gettime(CLOCK_MONOTONIC, &c_end0);
            if(signal)
            {   
                // clock_gettime(CLOCK_MONOTONIC, &c_start1);
                MPI_Recv(recfrnext.data(), next_size, MPI_FLOAT, rank+1, Multi_Data,
                        g_active_comm, MPI_STATUS_IGNORE);
                // clock_gettime(CLOCK_MONOTONIC, &c_end1);
                mergeSortedArrays(local_data, recfrnext, buffer);
                // clock_gettime(CLOCK_MONOTONIC, &c_start2);
                MPI_Send(&buffer[local_size], next_size, MPI_FLOAT, rank+1, Multi_Data, g_active_comm);
                // clock_gettime(CLOCK_MONOTONIC, &c_end2);
                for(int i=0; i<local_size; i++)
                    local_data[i] = buffer[i];
            }
        }
        else if(( (rank+count)%2 == 1) && rank != 0)
        {  
            // clock_gettime(CLOCK_MONOTONIC, &c_start3);
            MPI_Send(&leftmost_local, 1, MPI_FLOAT, rank-1, Single_Data, g_active_comm);
            MPI_Recv(&signal, 1, MPI_C_BOOL, rank-1, Bool, g_active_comm, MPI_STATUS_IGNORE);
            // clock_gettime(CLOCK_MONOTONIC, &c_end3);
            if(signal)
            {   
                // clock_gettime(CLOCK_MONOTONIC, &c_start4);
                MPI_Send(local_data.data(), local_size, MPI_FLOAT, rank-1, Multi_Data, g_active_comm);
                MPI_Recv(local_data.data(), local_size, MPI_FLOAT, rank-1, Multi_Data, g_active_comm, MPI_STATUS_IGNORE);
                // clock_gettime(CLOCK_MONOTONIC, &c_end4);
            }
        }
        
        //奇偶循环
        count = !count;
        //判断是否继续sort
        // clock_gettime(CLOCK_MONOTONIC, &c_start5);
        MPI_Reduce(&signal, &buff, 1, MPI_C_BOOL, MPI_LOR, 0, g_active_comm);
        MPI_Bcast(&buff, 1, MPI_C_BOOL, 0, g_active_comm);
        // clock_gettime(CLOCK_MONOTONIC, &c_end5);
        
        // temptimeCalculate(c_start0, c_end0, c_temp0);
        // temptimeCalculate(c_start3, c_end3, c_temp3);
        // temptimeCalculate(c_start5, c_end5, c_temp5);

        // if(signal)
        // {
        //     temptimeCalculate(c_start1, c_end1, c_temp1);
        //     temptimeCalculate(c_start2, c_end2, c_temp2);
        //     temptimeCalculate(c_start4, c_end4, c_temp4);
        //     c_time_used = c_temp0.tv_sec + c_temp1.tv_sec + c_temp2.tv_sec + c_temp3.tv_sec +
        //                 c_temp4.tv_sec + c_temp5.tv_sec + (double)(c_temp0.tv_nsec + c_temp1.tv_nsec+
        //                 c_temp2.tv_nsec + c_temp3.tv_nsec + c_temp4.tv_nsec + c_temp5.tv_nsec)/1e9;
        // }
        // else
        // {
        //     c_time_used = c_temp0.tv_sec + c_temp3.tv_sec + c_temp5.tv_sec +
        //                  (double)(c_temp0.tv_nsec + c_temp3.tv_nsec + c_temp5.tv_nsec)/1e9;
        // }
        // comm_time.push_back(c_time_used);

    }
}

void mergeSortedArrays(std::vector<float>& arr1, std::vector<float>& arr2, std::vector<float>& buffer) {
    int i = 0, j = 0;
    int arr1count = arr1.size();
    int arr2count = arr2.size();
    int count = 0;

    // 遍历两个数组并进行合并
    while (i < arr1count && j < arr2count) {
        if (arr1[i] < arr2[j]) {
            buffer[count++] = arr1[i];
            i++;
        } else {
            buffer[count++] = arr2[j];
            j++;
        }
    }
    // 将剩余元素添加到合并后的数组
    while (i < arr1count) {
        buffer[count++] = arr1[i];
        i++;
    }
    while (j < arr2count) {
        buffer[count++] = arr2[j];
        j++;
    }
}

// void temptimeCalculate(const struct timespec & start, const struct timespec & end, struct timespec & temp)
// {
//     if ((end.tv_nsec - start.tv_nsec) < 0) 
//     {
//        temp.tv_sec = end.tv_sec-start.tv_sec-1;
//        temp.tv_nsec = 1e9 + end.tv_nsec - start.tv_nsec;
//     } else {
//        temp.tv_sec = end.tv_sec - start.tv_sec;
//        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
//     }
// }

// void timeCalculate()
// {   
//     //total time calculate
//     temptimeCalculate(start, end, temp);
//     total_time_used = temp.tv_sec + (double)temp.tv_nsec/1e9;
    
//     //io time calculate
//     temptimeCalculate(i_start, i_end, i_temp);
//     temptimeCalculate(o_start, o_end, o_temp);
//     io_time_used = i_temp.tv_sec + o_temp.tv_sec + (double)(i_temp.tv_nsec + o_temp.tv_nsec)/1e9;

//     //comm time calculate
//     c_time_used = std::accumulate(comm_time.begin(), comm_time.end(), (double)0);
// }
