# Benchmarking Results

This section shows the benchmark runtimes and memory costs of our algorithms. 

Show the comparison with svd and rsvd based compression from tensor to mps and clearly see the advantage rsvd offers in terms of runtime at the cost of certain truncation errors. 

Show the benchmark runtime of the transform with and without preprocessing of your input data, for different signal types. The benchmark must run upto 60 output qubits in dt and zt, and run upto 30 qubits for qft. Compare the runtime of fftw also for the qft case (like chen and stoudenmire did). 

