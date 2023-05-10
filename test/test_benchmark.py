import numpy as np
import time

# Define the parameters
promoter_start_bin = 0
promoter_end_bin = 10000
enhancer_start_bin_index = 5000
enhancer_end_bin_index = 15000
re_pair_index = 0
pixel_index_array = [[] for i in range(100)]

# Benchmark the original code
start_time = time.time()
for promoter_bin_index in range(promoter_start_bin, promoter_end_bin + 1):
    for enhancer_bin_index in range(enhancer_start_bin_index, enhancer_end_bin_index + 1):
        pixel_index = (promoter_bin_index,enhancer_bin_index)
        pixel_index_array[re_pair_index].append(np.array([pixel_index[0],pixel_index[1]]))
end_time = time.time()
original_time = end_time - start_time

# Benchmark the NumPy vectorized code
start_time = time.time()
promoter_bin_indices = np.arange(promoter_start_bin, promoter_end_bin + 1)
enhancer_bin_indices = np.arange(enhancer_start_bin_index, enhancer_end_bin_index + 1)
pixel_indices = np.transpose([np.tile(promoter_bin_indices, len(enhancer_bin_indices)), np.repeat(enhancer_bin_indices, len(promoter_bin_indices))])
pixel_index_array[re_pair_index].append(pixel_indices)
end_time = time.time()
vectorized_time = end_time - start_time

print("Original code took", original_time, "seconds")
print("NumPy vectorized code took", vectorized_time, "seconds")
