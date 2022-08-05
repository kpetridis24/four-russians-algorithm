# Four-Russians-Algorithm
Boolean matrix multiplication accelerated by the four-Russians algorithm.

Boolean matrix multiplication has a default complexity of O(n^3), so using the Four-Russians-Algorithm
and multiplying by blocks of size t, we are able to drop the complexity to O(n^3/t). The idea is that
if we multiply Aik*Bkj the Cij will be a combination of some rows of B. These rows are indicated by the
corresponding subrow of A. For example lets say that we multiply the blocks:

          1 0  
     A1 = 1 1    and    B1 = 1 0 0 1 1 0
          0 1                0 0 1 0 0 1
     
The F-R algorithm takes advantage of the following fact. The rows of the C1 block, will be indicated by the 
A1's indexes. So as we see from the first row of A1, the C1's first row is the first row of B (1 0). The second
row is the "OR"(sum) of the two rows of B (1 1). And the C1's third row is the second row of B (0 1). So it is 
clear that we can pre-calculate all the possible row-sums of B, and looking at each row of A, choose the corresponding
row of our precomputed array. So the C1 block will be:

          1 0 0 1 1 0
     C1 = 1 0 1 1 1 1 
          0 0 1 0 0 1

Here i have to underline that such implementation would make sense in dense matrices.
