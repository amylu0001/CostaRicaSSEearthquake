To compile the file:
mpifort sub_comdun.f mod_dtrigreen.f90 m_calc_green.f90 calc_ds_new.f90 -o green_run

To calculate your own Green's function from a geometry file:
1. Prepare your geometry file, write it into a text file. The general format is:
Line 1: n_vertex, n_edge, n_cell   <--- total numbers of the vertex, edge, and cell
Line 2 to [1+n_vertex]: vertex(:,1), vertex(:,2), vertex(:,3) <--- coordinates of each vertex
Line [2+n_vertex] to [1+n_vertex+n_edge]: edge(:, 1), edge(:, 2) <--- edge, connection between two vertices
Line [2+n_vertex+n_edge] to [1+n_vertex+n_edge+n_cell]: cell(:, 1), cell(:, 2), cell(:, 3) <--- three vertice that forms a trigle

* In the edge and cell lines, you should write the number of the vertices. For example, vertices 1, 2, and 3 form a triangle. The vertices connect with each other to form 3 edges. And these 3 egdes forms a triangle. You should write:
Edge line:    1, 3
              1, 2
              2, 3
vertex line:  1, 2, 3


2. Enter the filename at line 14:
character(*),parameter        :: fname="costarica1kmv2.gts"

3. Enter the number of CPUs at line 22, make sure mod(n_cell,nproc)=0:
parameter (nproc=148)

4. Pick a stiffness matrix to compute at line 588:
To compute the matrix for the shear stress:
arr_out(j-(myid)*Nt,i) = - parm_miu/100 * dot_product(arr_cl_v2(:,3,j),matmul(sig33(:,:),arr_cl_v2(:,2,j)))
To compute the matrix for the normal stress:
arr_out(j-(myid)*Nt,i) = - parm_miu/100 * dot_product(arr_cl_v2(:,3,j),matmul(sig33(:,:),arr_cl_v2(:,3,j)))

The outputs are written into binary files:
Position file to record the center of each cell: position.bin
Stiffness matrix file: trigreen_0.bin trigreen_1.bin trigreen_2.bin ... trigreen_(num_CPU-1).bin
The number of files equals the number of CPUs

Before running the code for the sequence model, please put the shear stiffness files into a directory, name it, and write it on Line 3 in the "parameter.txt" file. Same for the normal stiffness files, write the directory name on Line 4.
For example:
./shear_stiffness/  
./normal_stiffness/
